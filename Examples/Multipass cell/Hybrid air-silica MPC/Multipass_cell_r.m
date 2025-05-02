% This code simulates (soliton) self-compression in a multipass cell.
%
%   Omar et al., "Hybrid air-bulk multi-pass cell compressor for high pulse energies with full spatio-temporal characterization," Opt. Express 32(8), 13235-13248 (2024)
%
% This script employs the radially-symmetric scheme of the UPPE code, 
% rather than a full x-y dimension.
%
% Current simulation bottleneck comes from the number of saved points per
% propagation, which limits the max step size of the adaptive-dz scheme.
% For testing, reduce "num_save" to improve the speed.

clearvars; close all;

addpath('../../../UPPE3D algorithm/','../../../user_helpers/'); % add package paths

%% Setup temporal/spectral parameters
Nt = 2^10; % the number of time points
time_window = 10; % ps
dt = time_window/Nt; % ps

%% Setup spatial parameters
% Information for the Hankel transform
Nr = 2^10; % the number of radial sampling points
r_max = 8e-3; % the maximum radius; half of the spatial window
kr_max = 2e4; % the maximum kr vector
FHATHA_energy_restoration = false; % don't enable constant operations of FHATHA's energy restoration function (see details of FHATHA in readme)

[r,kr,...
 dr,dkr,...
 l0,exp_prefactor,r2_prefactor,...
 ifftQ] = FHATHA_info(Nr,r_max,kr_max);

% Arrange required Hankel information into "sim" for radially-symmetric
% UPPE to use later.
sim.FHATHA = struct('r',r,'kr',kr,...
                    'dr',dr,'dkr',dkr,...
                    'l0',l0,'exp_prefactor',exp_prefactor,'r2_prefactor',r2_prefactor,...
                    'ifftQ',ifftQ,...
                    'energy_restoration',FHATHA_energy_restoration);

%% Setup simulation parameters
sim.lambda0 = 1037e-9; % m; the center wavelength
sim.gpuDevice.Index = 1; % choose which GPU to use if you have multiple GPUs: 1,2,3...

% Please check this function for details.
[fiber,sim] = load_default_UPPE3D_propagate([],sim); % load default parameters

% Setup general parameters
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Silica material properties
fiber_silica = fiber;
fiber_silica.material = 'silica'; % for finding the refractive index in UPPE3D_propagate()
% Sellmeier coefficients
[a,b] = Sellmeier_coefficients(fiber_silica.material);
% Calculate the index difference using the Sellmeier equation to generate n(lambda)
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));

% silica material properties
fiber_silica.n = n_from_Sellmeier(lambda/1e3); % refractive index
fiber_silica.n2 = 2.3e-20; % m^2/W; nonlinear refractive index

%% Configure gas parameters for the gas_info().
% The air properties (Raman, n, n2, etc.) is loaded in gas_info()
gas.temperature = 288; % K
gas.pressure = 1e5; % Pa
gas.material = {'air'};

% Load parameters based on the configured parameters
%
% gas.Ng - 1/m^3; gas number density
% gas.(gas.material).(Raman_type).(Raman_parameters)
[fiber_air,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

% Remove vibrational-Raman contributions
% Running with them requires a huge frequency window.
% Since they shouldn't play out significantly here in the pulse-compression
% demonstration, we turn it off so that we can use a small frequency window
% with fewer sampling points.
gas.N2.V.preR = gas.N2.V.preR*0;
gas.O2.V.preR = gas.O2.V.preR*0;

%% Setup MPC parameters
num_roundtrip = 23;
MPC_length = 1.34; % m
focal_R = 2; % m; radius of curvature of the mirror
focal_length = focal_R/2; % m
silica_thickness = 1e-3; % m
silica_to_mirror = 0.3; % m

plate_z = MPC_length/2;
for i = 1:num_roundtrip*6
    switch mod(i,6)
        case {0,1}
            zi = MPC_length - silica_to_mirror - silica_thickness;
        case {2,5}
            zi = silica_thickness;
        case {3,4}
            zi = silica_to_mirror;
    end
    plate_z = [plate_z;plate_z(end)+zi];
end

%% Initial condition
MFD0 = 1.4e-3; % m; mode-field diameter
tfwhm = 0.220; % ps
energy = 400e3; % nJ
initial_condition = build_3Dgaussian_r(MFD0, tfwhm, time_window, energy, Nt, r);

% In the paper, they chirped the pulse to 276 fs.
func = calc_chirp;
chirp_sign = 1; % for positive chirp (normal dispersive delay)
duration = 0.276; % ps
[chirp,chirped_pulse] = func.General( duration,ifftshift(2*pi*f,1),ifft(initial_condition.field,[],1),chirp_sign );
initial_condition.field = chirped_pulse;

%%
% =========================================================================
% Start the simulation of MPC
% =========================================================================
num_save = 5;
z_all = zeros(num_save,6*num_roundtrip+2);
MFD_all = zeros(num_save,6*num_roundtrip+2);

Frame(num_save,6*num_roundtrip+1) = struct('cdata',[],'colormap',[]);

%% Before hitting the first mirror
fiber_air.L0 = MPC_length/2; % m

sim.save_period = fiber_air.L0/num_save;

% Show initial k space
A0_H = 2*pi*FHATHA(squeeze(initial_condition.field(floor(Nt/2)+1,:,end)),...
                   r_max,...
                   r,kr,...
                   dr,dkr,...
                   l0,exp_prefactor,r2_prefactor,...
                   ifftQ);
fig_k = figure;
plot(kr/1e3,abs(A0_H).^2,'linewidth',2,'Color','r');
xlabel('k_r (2\pi/mm)');
set(gca,'fontsize',20);
title('Initial k space');
% Plot the 2D field with pcolor
fig_k2 = radialPcolor(kr/1e3,abs(A0_H).^2);
xlabel('k_x (2\pi/mm)');
ylabel('k_y (2\pi/mm)');
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Initial k space');

% Simulate the propagation
prop_output = UPPE3D_propagate(fiber_air,initial_condition,sim,gas);

% Calculate MFD during propagation
z_all(:,1) = prop_output.z(2:end); % the first one is the input z=0, so ignore it
remove_noise_model = 2; % see calcMFD for details
MFD = calcMFD_r(squeeze(sqrt(sum(abs(prop_output.field).^2,1))),r,remove_noise_model)'*1e3; % mm
MFD_all(:,1) = MFD(2:end); % the first one is the input, so ignore it
fprintf('air %u\n',0);
num_to_plot = num_save;
fig = plotter_r([],...
                prop_output.field,...
                z_all(1:num_to_plot),MFD_all(1:num_to_plot),...
                num_save,...
                r,lambda);

% Animation
Frame(:,1) = animator_r(Frame(:,1),...
                        prop_output.field,...
                        z_all(1:num_to_plot),MFD_all(1:num_to_plot),0,...
                        num_save,...
                        Nt,dt,r,lambda,...
                        plate_z);

%% MPC
for i = 1+(1:num_roundtrip*6)
    initial_condition.field = prop_output.field(:,:,:,end);

    if ismember(mod(i-1,6),[1,4])
        Ef_reflected = add_thin_lens_phase_r(ifft(initial_condition.field,[],1),initial_condition.r,ifftshift(lambda,1)/1e9,focal_length);
        initial_condition.field = fft(Ef_reflected,[],1);
    end
    
    switch mod(i-1,6)
        case {0,1}
            fiber_air.L0 = MPC_length - silica_to_mirror - silica_thickness;
        case {2,5}
            fiber_silica.L0 = silica_thickness;
        case {3,4}
            fiber_air.L0 = silica_to_mirror;
    end

    % Show initial k space
    close(fig_k,fig_k2);
    A0_H = 2*pi*FHATHA(squeeze(initial_condition.field(floor(Nt/2)+1,:,end)),...
                       r_max,...
                       r,kr,...
                       dr,dkr,...
                       l0,exp_prefactor,r2_prefactor,...
                       ifftQ);
    fig_k = figure;
    plot(kr/1e3,abs(A0_H).^2,'linewidth',2,'Color','r');
    xlabel('k_r (2\pi/mm)');
    set(gca,'fontsize',20);
    title('Initial k space');
    % Plot the 2D field with pcolor
    fig_k2 = radialPcolor(kr/1e3,abs(A0_H).^2);
    xlabel('k_x (2\pi/mm)');
    ylabel('k_y (2\pi/mm)');
    set(gca,'fontsize',20);
    daspect([1 1 1]); % make aspect ratio = 1
    title('Initial k space');
    
    % Simulate the propagation
    switch mod(i-1,6)
        case {0,1,3,4} % in air
            sim.save_period = fiber_air.L0/num_save;
            prop_output = UPPE3D_propagate(fiber_air,initial_condition,sim,gas);
        case {2,5} % in silica
            sim.save_period = fiber_silica.L0/num_save;
            prop_output = UPPE3D_propagate(fiber_silica,initial_condition,sim);
    end
    
    % Calculate MFD during propagation
    z_all(:,i) = prop_output.z(2:end) + z_all(end,i-1); % the first one is the input z=0, so ignore it
    MFD = calcMFD_r(squeeze(sqrt(sum(abs(prop_output.field).^2,1))),r,remove_noise_model)'*1e3; % mm
    MFD_all(:,i) = MFD(2:end); % the first one is the input, so ignore it
    switch mod(i-1,6)
        case {0,1,3,4}
            fprintf('air %u\n',floor((i-1)/6));
        case {2,5}
            fprintf('silica %u\n',floor((i-1)/6));
    end
    num_to_plot = i*num_save;
    fig = plotter_r(fig,... ...
                    prop_output.field,...
                    z_all(1:num_to_plot),MFD_all(1:num_to_plot),...
                    num_save,...
                    r,lambda);

    % Animation
    Frame(:,i) = animator_r(Frame(:,i),...
                            prop_output.field,...
                            z_all(1:num_to_plot),MFD_all(1:num_to_plot),(i-1)*num_save,...
                            num_save,...
                            Nt,dt,r,lambda,...
                            plate_z);
end

% Movie
implay(Frame(:),num_save);
exportVideo = VideoWriter('MPC_r');
exportVideo.FrameRate = num_save;
open(exportVideo);
writeVideo(exportVideo, Frame(:));
close(exportVideo);

%% Dechirping
% Use the center field as the whole temporal profile for dechirping. This
% ignores spatiotemporal coupling.
output_mode_field = prop_output.field(:,1,end)*sqrt(pi*(MFD(end)/2*1e-3)^2);

theta_in = pi/6;
wavelength0 = sim.lambda0*1e9;
grating_spacing = 1e-6;
pulse_compressor( 'Treacy-t',theta_in,wavelength0,t,output_mode_field,grating_spacing,true );