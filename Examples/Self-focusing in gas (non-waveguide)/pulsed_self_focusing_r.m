% This code simulates self-focusing of a Gaussian pulse in N2.
%
% Gas-filled 3D-UPPE employs the radially-symmetric scheme of the UPPE code.

close all; clearvars;

addpath('../../UPPE3D algorithm/','../../user_helpers/');

%% Setup temporal/spectral parameters
Nt = 2^10; % the number of time points
time_window = 20; % ps
dt = time_window/Nt; % ps

%% Setup spatial parameters
% Information for the Hankel transform
Nr = 2^11; % the number of radial sampling points
r_max = 5e-4; % the maximum radius; half of the spatial window
kr_max = 5e5; % the maximum kr vector
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

%% Setup fiber parameters
sim.lambda0 = 1030e-9; % the center wavelength
sim.gpuDevice.Index = 1;

fiber.L0 = 0.002;
num_save = 30;
sim.save_period = fiber.L0/num_save;

% Load default parameters like
%
% load fiber.betas and fiber.SR based on your multimode folder above
% sim.f0 = 3e5/sim.lambda0; THz
% sim.save_period = 0; Save only the field at the beginning and the end fiber
% sim.scalar = true; Use scalar propagation
% sim.adaptive_dz.threshold = 1e-5; the threshold of the adaptive-step method
% sim.gpu_yes = true; Use GPU
% sim.include_Raman = 1; Consider the Raman
% sim.pulse_centering = true; Always shift the pulse to the center of the time window
% sim.gpuDevice.Index = 1; Use the GPU device 1
% sim.progress_bar = true; Show the progress bar
% sim.progress_bar_name = ''; Empty name for the progress bar
% sim.cuda_dir_path = 'UPPE3D/cuda'; Where the cuda files are
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_UPPE3D_propagate(fiber,sim); % load default parameters

% Setup general parameters
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Configure gas parameters for the gas_info().
gas.temperature = 288; % K
gas.pressure = 10e5; % Pa
gas.material = {'N2'};

sim.photoionization_model = true; % enable photoionization

% Load parameters based on the configured parameters
%
% gas.Ng - 1/m^3; gas number density
% gas.(gas.material).(Raman_type).(Raman_parameters)
[fiber,sim,gas] = gas_info(fiber,sim,gas,lambda*1e-9);

% Remove vibrational-Raman contributions
% Running with them requires a huge frequency window.
% Since they shouldn't play out significantly here in the pulse-compression
% demonstration, we turn it off so that we can use a small frequency window
% with fewer sampling points.
gas.N2.V.preR = gas.N2.V.preR*0;

%% Initial condition
MFD0 = 100e-6; % m
tfwhm = 1; % ps
energy = 3e6; % nJ
initial_condition = build_3Dgaussian_r(MFD0, tfwhm, time_window, energy, Nt, r);

%% Show initial spaces
% Show initial real space
figure;
plot(r*1e6,abs(squeeze(initial_condition.field(floor(Nt/2)+1,:))).^2,'linewidth',2,'Color','b');
xlabel('r (\mum)');
xlim([0,100]);
set(gca,'fontsize',20);
title('Initial real space');
% Plot the 2D field with pcolor
radialPcolor(r*1e6,abs(squeeze(initial_condition.field(floor(Nt/2)+1,:,end))).^2);
xlabel('x (\mum)');
ylabel('y (\mum)');
xlim([-100,100]);
ylim([-100,100]);
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Initial real space');

A0_H = 2*pi*FHATHA(squeeze(initial_condition.field(floor(Nt/2)+1,:)),...
                   r_max,...
                   r,kr,...
                   dr,dkr,...
                   l0,exp_prefactor,r2_prefactor,...
                   ifftQ);

% Show initial k space
figure;
plot(kr/1e6,abs(A0_H).^2,'linewidth',2,'Color','r');
xlabel('k_r (2\pi/\mum)');
set(gca,'fontsize',20);
title('Initial k space');
% Plot the 2D field with pcolor
radialPcolor(kr/1e6,abs(A0_H).^2);
xlabel('k_x (2\pi/\mum)');
ylabel('k_y (2\pi/\mum)');
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Initial k space');

%% Propagate
prop_output = UPPE3D_propagate(fiber,initial_condition,sim,gas);

%% Results
MFD = calcMFD_r(squeeze(prop_output.field(floor(Nt/2)+1,:,:)),r)'*1e6; % um
optical_energy = squeeze(sum(2*pi*trapz(r,abs(prop_output.field).^2.*r,2),1)*dt/1e3); % nJ

%% Plot
% Show final real space
figure;
plot(r*1e6,abs(squeeze(prop_output.field(floor(Nt/2)+1,:,end))).^2,'linewidth',2,'Color','b');
xlabel('r (\mum)');
xlim([0,100]);
set(gca,'fontsize',20);
title('Final real space');
% Plot the 2D field with pcolor
radialPcolor(r*1e6,abs(squeeze(prop_output.field(floor(Nt/2)+1,:,end))).^2);
xlabel('x (\mum)');
ylabel('y (\mum)');
xlim([-100,100]);
ylim([-100,100]);
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Final real space');
% Show final k space
A_H = 2*pi*FHATHA(squeeze(prop_output.field(floor(Nt/2)+1,:,end)),...
                  r_max,...
                  r,kr,...
                  dr,dkr,...
                  l0,exp_prefactor,r2_prefactor,...
                  ifftQ);
figure;
plot(kr/1e6,abs(A_H).^2,'linewidth',2,'Color','r');
xlabel('k_r (2\pi/\mum)');
set(gca,'fontsize',20);
title('Final k space');
% Plot the 2D field with pcolor
radialPcolor(kr/1e6,abs(A_H).^2);
xlabel('k_x (2\pi/\mum)');
ylabel('k_y (2\pi/\mum)');
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Final k space');

% Plot MFD
figure;
plot(prop_output.z*1e2,MFD,'linewidth',2,'Color','k');
xlabel('Propagation distance (cm)');
ylabel('MFD (\mum)');
set(gca,'fontsize',20);

% Energy
figure;
plot(prop_output.z,optical_energy,'linewidth',2,'Color','b');
xlabel('Propagation distance (m)');
ylabel('Power (nJ)');
set(gca,'fontsize',20);

% Plot spatial evolution
plot_spatial_evolution_r(r*1e6,prop_output.z,squeeze(abs(prop_output.field(floor(Nt/2)+1,:,:)).^2));
ylabel('X (\mum)');
ylim([-100,100]);

% Movie
Frame = pulsed_animator_r(prop_output.field,....
                          prop_output.z,MFD,...
                          Nt,dt,r,lambda,...
                          fiber.L0);
implay(Frame(:),20);
exportVideo = VideoWriter('pulsed_self-focusing');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame(:));
close(exportVideo);

%save('self_focusing_r.mat');