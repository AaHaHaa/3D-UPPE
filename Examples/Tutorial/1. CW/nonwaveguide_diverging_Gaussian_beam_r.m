% This code simulates a diverging Gaussian beam in air. The simulated MFD
% is compared to the theoretical values.
%
% This script employs the radially-symmetric scheme of the UPPE code, 
% rather than a full x-y dimension.
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
% Information for the Hankel transform
Nr = 2^8; % the number of radial sampling points
r_max = 4e-3; % the maximum radius; half of the spatial window
kr_max = 3e4; % the maximum kr vector
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
fiber.n = 1; % air
fiber.n2 = 0; % no nonlinearity

fiber.L0 = 1;
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

% General parameters
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Initial condition
% Pulse duration and energy itself have no meaning here. Only the peak
% power is retained as the CW power in building the initial profile.
% However, still ensure that pulse duration is around 5-10x smaller than
% the time window for correct generation of the CW profile.
MFD0 = 500e-6; % m
tfwhm = 1; % ps; determine the peak power, which is average power in CW
energy = 1e-3; % nJ; determine the peak power, which is average power in CW
initial_condition = build_3Dgaussian_r(MFD0, tfwhm, time_window, energy, Nt, r);

%% Show initial spaces
% Show initial real space
figure;
plot(r*1e3,abs(squeeze(initial_condition.field(floor(Nt/2)+1,:))).^2,'linewidth',2,'Color','b');
xlabel('r (mm)');
set(gca,'fontsize',20);
title('Initial real space');
% Plot the 2D field with pcolor
radialPcolor(r*1e3,abs(squeeze(initial_condition.field(floor(Nt/2)+1,:,end))).^2);
xlabel('r (mm)');
ylabel('r (mm)');
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
plot(kr/1e3,abs(A0_H).^2,'linewidth',2,'Color','r');
xlabel('k_r (2\pi/mm)');
set(gca,'fontsize',20);
title('Initial k space');
% Plot the 2D field with pcolor
radialPcolor(kr/1e3,abs(A0_H).^2);
xlabel('k_x (2\pi/mm)');
ylabel('k_y (2\pi/mm)');
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Initial k space');

%% Propagate
prop_output = UPPE3D_propagate(fiber,initial_condition,sim);

%% Results
MFD = calcMFD_r(squeeze(prop_output.field(floor(Nt/2)+1,:,:)),r)'*1e3; % mm
optical_power = squeeze(sum(2*pi*trapz(r,abs(prop_output.field).^2.*r,2),1)); % W

%% Theoretical Gaussian propagation
w0 = MFD0/2;
zR = pi*w0^2/sim.lambda0; % Raylength length
MFD_theory = MFD0*sqrt(1+(squeeze(prop_output.z)/zR).^2)*1e3; % mm

%% Plot
% Show final real space
figure;
plot(r*1e3,abs(squeeze(prop_output.field(floor(Nt/2)+1,:,end))).^2,'linewidth',2,'Color','b');
xlabel('r (mm)');
set(gca,'fontsize',20);
title('Final real space');
% Plot the 2D field with pcolor
radialPcolor(r*1e3,abs(squeeze(prop_output.field(floor(Nt/2)+1,:,end))).^2);
xlabel('r (mm)');
ylabel('r (mm)');
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
plot(kr/1e3,abs(A_H).^2,'linewidth',2,'Color','r');
xlabel('k_r (2\pi/mm)');
set(gca,'fontsize',20);
title('Final k space');
% Plot the 2D field with pcolor
radialPcolor(kr/1e3,abs(A_H).^2);
xlabel('k_x (2\pi/mm)');
ylabel('k_y (2\pi/mm)');
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
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