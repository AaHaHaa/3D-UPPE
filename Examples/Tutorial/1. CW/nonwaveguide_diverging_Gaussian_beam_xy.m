% This code simulates a diverging Gaussian beam in air. The simulated MFD
% is compared to the theoretical values.
%
% This script employs the 3D-UPPE that uses full x-y dimension. For
% more-efficient modeling, pelase see its radially-symmetric version.
%
% In addition, this script simulates with continuous-wave (CW) light, which
% focuses only on the evolution of its spatial profile through the
% diffraction effect.
%
% To run with CW, make Nt = 1.

close all; clearvars;

addpath('../../../UPPE3D algorithm/','../../../user_helpers/');

%% Setup temporal/spectral parameters
Nt = 1; % the number of time points
time_window = 10; % ps; not important here due to CW
dt = time_window/Nt; % ps; not important here due to CW

%% Setup spatial parameters
Nx = 2^7; % the number of spatial points
spatial_window = 10e-3; % m
dx = spatial_window/Nx;

%% Setup fiber parameters
sim.lambda0 = 1030e-9; % the center wavelength
fiber.n = ones(1,Nx,Nx); % air
fiber.n2 = 0; % no nonlinearity

fiber.L0 = 1; % m
num_save = 100;
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
f = sim.f0+(-floor(Nt/2):floor((Nt-1)/2))'/(Nt*dt); % THz
t = (-floor(Nt/2):floor((Nt-1)/2))'*dt; % ps; it is meaningless in this code since Nt = 1 (CW case)
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

x = (-Nx/2:Nx/2-1)*dx*1e3; % mm
kx = 2*pi*(-Nx/2:Nx/2-1)/spatial_window*1e-3; % 2*pi/mm

%% Initial condition
% Pulse duration and energy itself have no meaning here. Only the peak
% power is retained as the CW power in building the initial profile.
% However, still ensure that pulse duration is around 5-10x smaller than
% the time window for correct generation of the CW profile.
MFD0 = 500e-6; % m
tfwhm = 1; % ps; determine the peak power, which is average power in CW
energy = 1e-3; % nJ; determine the peak power, which is average power in CW
initial_condition = build_3Dgaussian_xy(MFD0, spatial_window, tfwhm, time_window, energy, Nt, Nx);

%% Show initial spaces
% Show initial real space
figure;
pcolor(x,x,abs(squeeze(initial_condition.field(floor(Nt/2)+1,:,:))).^2); colormap(jet);
shading interp;colormap(jet);colorbar;
xlabel('x (mm)');
ylabel('y (mm)');
daspect([1 1 1]); % make aspect ratio = 1
set(gca,'fontsize',20);
title('Initial real space');
% Show initial k space
figure;
pcolor(kx,kx,abs(fftshift(fft(fft(squeeze(initial_condition.field(floor(Nt/2)+1,:,:)),[],1),[],2))).^2); colormap(jet);
shading interp;colormap(jet);colorbar;
xlabel('k_x (2\pi/mm)');
ylabel('k_y (2\pi/mm)');
daspect([1 1 1]); % make aspect ratio = 1
set(gca,'fontsize',20);
title('Initial k space');

%% Propagate
prop_output = UPPE3D_propagate(fiber,initial_condition,sim);

%% Results
MFD = squeeze(calcMFD_xy(squeeze(prop_output.field(floor(Nt/2)+1,:,:,:)),spatial_window))*1e3; % mm
optical_power = squeeze(sum(abs(prop_output.field).^2,[1,2,3]))*dx^2; % W

%% Theoretical Gaussian propagation
w0 = MFD0/2;
zR = pi*w0^2/sim.lambda0; % Raylength length
MFD_theory = MFD0*sqrt(1+(squeeze(prop_output.z)/zR).^2)*1e3; % mm

%% Plot
% Show final real space
figure;
pcolor(x,x,abs(squeeze(prop_output.field(floor(Nt/2)+1,:,:,end))).^2); colormap(jet);
shading interp;colormap(jet);colorbar;
xlabel('x (mm)');
ylabel('y (mm)');
daspect([1 1 1]); % make aspect ratio = 1
set(gca,'fontsize',20);
title('Final real space');
% Show final k space
figure;
pcolor(kx,kx,abs(fftshift(fft(fft(squeeze(prop_output.field(floor(Nt/2)+1,:,:,end)),[],1),[],2))).^2); colormap(jet);
shading interp;colormap(jet);colorbar;
xlabel('x (2\pi/mm)');
ylabel('y (2\pi/mm)');
daspect([1 1 1]); % make aspect ratio = 1
set(gca,'fontsize',20);
title('Final k space');

% Plot MFD
figure;
h = plot(prop_output.z,[MFD,MFD_theory],'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r');
set(gca,'fontsize',20);
xlabel('Propagation distance (m)');
ylabel('MFD (mm)');
l = legend('Simulated','Calculated'); set(l,'location','northwest');

% Power
% Check that whether it's conserved
figure;
plot(prop_output.z,optical_power,'linewidth',2,'Color','b');
xlabel('Propagation distance (m)');
ylabel('Power (W)');
set(gca,'fontsize',20);