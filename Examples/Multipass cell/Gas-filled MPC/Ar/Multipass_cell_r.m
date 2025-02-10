% This code simulates a multipass cell and aims to duplicate Figs. 1 and 2 of the reference paper:
%
%   Ueffing et al., "Nonlinear pulse compression in a gas-filled multipass cell," Opt. Lett. 43(9), 2070-2073 (2018)
%
% The generated 35-fs dechirped pulse matches well with the experimental 37 fs in the paper.
%
% Gas-filled 3D-UPPE employs the radially-symmetric scheme of the UPPE code.

clearvars; close all;

addpath('../../../../UPPE3D algorithm/','../../../../user_helpers/'); % add package paths

%% Setup temporal/spectral parameters
Nt = 2^8; % the number of time points
time_window = 1; % ps
dt = time_window/Nt; % ps

%% Setup spatial parameters
% Information for the Hankel transform
Nr = 2^11; % the number of radial sampling points
r_max = 1e-2; % the maximum radius; half of the spatial window
kr_max = 5e4; % the maximum kr vector
FHATHA_energy_restoration = false; % don't enable constant operations of FHATHA's energy restoration function (see details of FHATHA in readme)

[r,kr,...
 dr,dkr,...
 l0,exp_prefactor,n2_prefactor,...
 ifftQ] = FHATHA_info(Nr,r_max,kr_max);

% Arrange required Hankel information into "sim" for radially-symmetric
% UPPE to use later.
sim.FHATHA = struct('r',r,'kr',kr,...
                    'dr',dr,'dkr',dkr,...
                    'l0',l0,'exp_prefactor',exp_prefactor,'n2_prefactor',n2_prefactor,...
                    'ifftQ',ifftQ,...
                    'energy_restoration',FHATHA_energy_restoration);

%% Setup simulation parameters
sim.lambda0 = 1030e-9; % m; the center wavelength
%sim.photoionization_model = true;
sim.gpuDevice.Index = 1; % choose which GPU to use if you have multiple GPUs: 1,2,3...

% Please check this function for details.
[fiber,sim] = load_default_UPPE3D_propagate([],sim); % load default parameters

% Setup general parameters
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Configure gas parameters for the gas_info().
gas.temperature = 288; % K
gas.pressure = 0.265*1e5; % Pa
gas.material = 'Ar';

% Load parameters based on the configured parameters
%
% gas.Ng - 1/m^3; gas number density
% gas.(gas.material).(Raman_type).(Raman_parameters)
[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

%% Setup MPC parameters
num_roundtrip = 24;
MPC_length = 1; % m
focal_R = 1; % m; radius of curvature of the mirror
focal_length = focal_R/2; % m

plate_z = MPC_length/2 + (1:num_roundtrip*2)*MPC_length;

%% Initial condition
MFD0 = 3.2e-3;%2*((focal_R^2-(MPC_length-focal_R)^2)*sim.lambda0^2/4/pi^2)^(1/4); % from Eq. (4) in the reference paper
tfwhm = 0.21; % ps
energy = 2e6; % nJ
initial_condition = build_3Dgaussian_r(MFD0, tfwhm, time_window, energy, Nt, r);

%%
% =========================================================================
% Start the simulation of MPC
% =========================================================================
num_save = 10;
z_all = zeros(num_save,2*num_roundtrip+1);
MFD_all = zeros(num_save,2*num_roundtrip+1);

Frame(num_save,2*num_roundtrip+1) = struct('cdata',[],'colormap',[]);

%% Before hitting the first mirror
fiber.L0 = MPC_length; % m
sim.save_period = fiber.L0/num_save;

Ef_reflected = add_thin_lens_phase_r(ifft(initial_condition.field,[],1),initial_condition.r,ifftshift(lambda,1)/1e9,focal_length);
initial_condition.field = fft(Ef_reflected,[],1);

% Show initial k space
A0_H = 2*pi*FHATHA(squeeze(initial_condition.field(floor(Nt/2)+1,:,end)),...
                   r_max,...
                   r,kr,...
                   dr,dkr,...
                   l0,exp_prefactor,n2_prefactor,...
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
prop_output = UPPE3D_propagate(fiber,initial_condition,sim,gas);

% Calculate MFD during propagation
z_all(:,1) = prop_output.z(2:end); % the first one is the input z=0, so ignore it
remove_noise_model = 2; % see calcMFD for details
MFD = calcMFD_r(squeeze(sqrt(sum(abs(prop_output.field).^2,1))),r,remove_noise_model)'*1e3; % mm
MFD_all(:,1) = MFD(2:end); % the first one is the input, so ignore it
fprintf('Ar %u\n',0);
num_to_plot = num_save;
fig = plotter_r([],...
                prop_output.field,...
                z_all(1:num_to_plot),MFD_all(1:num_to_plot),...
                r,lambda);

% Animation
Frame(:,1) = animator_r(Frame(:,1),...
                        prop_output.field,...
                        z_all(1:num_to_plot),MFD_all(1:num_to_plot),0,...
                        Nt,dt,r,lambda,...
                        plate_z);

%% MPC
fiber.L0 = MPC_length;
for i = 1+(1:num_roundtrip*2)
    initial_condition.field = prop_output.field(:,:,:,end);

    Ef_reflected = add_thin_lens_phase_r(ifft(initial_condition.field,[],1),initial_condition.r,ifftshift(lambda,1)/1e9,focal_length);
    initial_condition.field = fft(Ef_reflected,[],1);
    
    sim.save_period = fiber.L0/num_save;

    % Show initial k space
    close(fig_k,fig_k2);
    A0_H = 2*pi*FHATHA(squeeze(initial_condition.field(floor(Nt/2)+1,:,end)),...
                       r_max,...
                       r,kr,...
                       dr,dkr,...
                       l0,exp_prefactor,n2_prefactor,...
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
    prop_output = UPPE3D_propagate(fiber,initial_condition,sim,gas);
    
    % Calculate MFD during propagation
    z_all(:,i) = prop_output.z(2:end) + z_all(end,i-1); % the first one is the input z=0, so ignore it
    MFD = calcMFD_r(squeeze(sqrt(sum(abs(prop_output.field).^2,1))),r,remove_noise_model)'*1e3; % mm
    MFD_all(:,i) = MFD(2:end); % the first one is the input, so ignore it
    
    fprintf('Ar %u\n',floor((i-1)/2));

    num_to_plot = i*num_save;
    fig = plotter_r(fig,... ...
                    prop_output.field,...
                    z_all(1:num_to_plot),MFD_all(1:num_to_plot),...
                    r,lambda);

    % Animation
    Frame(:,i) = animator_r(Frame(:,i),...
                            prop_output.field,...
                            z_all(1:num_to_plot),MFD_all(1:num_to_plot),(i-1)*num_save,...
                            Nt,dt,r,lambda,...
                            plate_z);
end

% Movie
implay(Frame(:),30);
exportVideo = VideoWriter('MPC_r');
exportVideo.FrameRate = 30;
open(exportVideo);
writeVideo(exportVideo, Frame(:));
close(exportVideo);

%% Dechirping
% Use the center field as the whole temporal profile for dechirping. This
% ignores spatiotemporal coupling.
output_mode_field = prop_output.field(:,1,end)*sqrt(pi*(MFD(end)/2*1e-3)^2);

theta_in = pi/6;
wavelength0 = sim.lambda0*1e6;
grating_spacing = 1e-6;
pulse_compressor( 'Treacy-t',theta_in,wavelength0,t,output_mode_field,grating_spacing,true );