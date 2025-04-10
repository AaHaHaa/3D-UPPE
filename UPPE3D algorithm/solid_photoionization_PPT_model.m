function [ne,DneDt] = solid_photoionization_PPT_model(Et, ionization_energy, f0, dt, Ne, n_pulse,...
                                                      me, tau_c, tau_r,...
                                                      S_Q, ellip_x, ellipK, ellipE, erfi_x, erfi_y,...
                                                      ellipticity,...
                                                      F_op)
%SOLID_PHOTOIONIZATION_PPT_MODEL This code computes the ionized free electrons
%with the Perelomov-Popov-Terent'ev (PPT) photoionization model. It's used
%as a helper function in the gas-UPPE model.
%   
% Input:
%
%   Et: (Nt,...); the electric field in the time domain
%   ionization_energy: scalar (J)
%   f0: center frequency of the frequency window
%   dt: scalar (ps)
%   Ne: the number density of neutral electrons (1/m^3)
%   n_pulse: temporally-averaged refractive index
%   me: reduced electron mass (kg)
%   tau_c: collision time (s)
%   tau_r: recombination time (s)
%   n_Q: the number of terms to be summed in ionization rate equation
%   ellip_x, ellipK, ellipE: (Nt,1); the lookup table for the complete elliptic integral of the first and second kinds
%   erfi_x, erfi_y: (Nt,1); the lookup table for the imaginary error function, erfi()
%   ellipticity: a scalar; only ellipticity=0 is allowed
%   F_op: Fourier-transform operators
%
% Output:
%   ne: (Nt,1); the free electron number density (1/m^3)
%   DneDt: (Nt,1); the time derivative of ne (1/m^3/s)
%
% Reference:
%   Couairon et al., "Filamentation and damage in fused silica induced by
%   tightly focused femtosecond laser pulses," Phys. Rev. B 71, 125435
%   (2005)
% Note that the paper by Couairon et al. has an error. The ionization rate
% W should be multiplied by (1-ne/Ng).

% The polarization dimension is in the 4th in radially-symmetric scheme and
% in 5th dimension in the xy scheme. Here, I kinda try to get them when I
% don't know what exactly the scheme the code is using.
[Nt,Npr,Npxy] = size(Et,[1,4,5]);
if Npr ~= 1 || Npxy ~= 1 || ellipticity ~= 0
    error('Photoionization model works only for linearly polarized fields.');
end

% Find instantaneous frequency of the pulse
pulse_phase = unwrap(angle(Et));
% Smoothing is required; otherwise, it's too noisy such that an erroneous
% high-frequency number is calculated
smooth_kernel = ones(floor(Nt/100),1)/floor(Nt/100);
Nk = length(smooth_kernel);
conv_length = Nk + Nt - 1;
smooth_kernel = F_op.Ff(smooth_kernel,conv_length);
pulse_phase = F_op.Ff(pulse_phase,conv_length);
pulse_phase = F_op.iFf(pulse_phase.*smooth_kernel,[]);
pulse_phase = real(pulse_phase(floor(Nk/2) : end-floor(Nk/2),:,:));
% conv() works for only a single-column vector
%pulse_phase = conv(pulse_phase,ones(floor(Nt/100),1)/floor(Nt/100),'same'); 
omega_pulse = -(pulse_phase(3:end,:,:)-pulse_phase(1:end-2,:,:))/(2*dt)+2*pi*f0; % THz; I use "central difference" to calculate the slope here
omega_pulse = cat(1,omega_pulse(1,:,:),omega_pulse,omega_pulse(end,:,:))*1e12; % Hz

e = 1.60217663e-19; % Coulomb
permittivity0 = 8.85418782e-12; % m^(-3)/kg*s^4*A^2
c = 299792458; % m/s
h = 6.62607015e-34; % m^2*kg/s
hbar = h/2/pi;

I = abs(Et).^2; % intensity; W/m^2

% This modification is for I too small to avoid spurious result in 1./I computation
I(I<max(I(:))/1e5) = max(I(:))/1e5;

ponderomotive_energy = e^2/2/me/permittivity0/c*I./omega_pulse.^2; % unlike gases, solid's refractive index isn't in Keldysh parameter
Keldysh_parameter = sqrt(ionization_energy/2./ponderomotive_energy);

Xi = 1./(1 + (Keldysh_parameter/sqrt(2)).^2); Xi(Xi>ellip_x(end)) = ellip_x(end);
Gamma = (Keldysh_parameter/sqrt(2)).^2.*Xi; Gamma(Gamma>ellip_x(end)) = ellip_x(end);
ellipKGamma = interp1(ellip_x,ellipK,Gamma);
ellipEGamma = interp1(ellip_x,ellipE,Gamma);
ellipKXi = interp1(ellip_x,ellipK,Xi);
ellipEXi = interp1(ellip_x,ellipE,Xi);
x = 2/pi*ionization_energy/hbar./omega_pulse.*ellipEXi./sqrt(Gamma);
v = ceil(x) - x;
alpha = pi*(ellipKGamma-ellipEGamma)./ellipEXi;
beta = pi^2/4./(ellipKXi.*ellipEXi);
Q = 0;
for S = 0:S_Q
    erfix = interp1(erfi_x,erfi_y,sqrt(beta.*(S+2*v)));
    Q = Q + exp(-S*(alpha+beta)).*erfix;
end
Q = pi/2/sqrt(2)./sqrt(ellipKXi).*exp(-2*beta.*v).*Q;

% the Keldysh ionization rate
W = 2*omega_pulse/9/pi.*(omega_pulse.*me/hbar./sqrt(Gamma)).^(3/2).*Q.*exp(-alpha.*ceil(x));

W(I<max(I(:))/1e4) = 0; % avoid non-negligible W when there is no or weak field due to numerical precision error

wtc = omega_pulse*tau_c;
sigma_IB = e^2./n_pulse./omega_pulse/c/permittivity0/me.*wtc./(1+wtc.^2);
P = cumsum((-W/Ne + sigma_IB/ionization_energy.*I - 1/tau_r)*(dt*1e-12),1);
ne = exp(P).*cumsum(W.*exp(-P),1)*(dt*1e-12);
if max(ne(:))/Ne > 0.3
    error('photoionization_PPT_model:neError',...
          ['Ion density is too high, reaching %.3f%% of total gas density.\n',...
           'Current implementation supports only perturbative ionization.\n',...
           'Please reduce the peak power or anything that reduces the ionization effect.'],max(ne(:))/Ne*100);
end

DneDt = W.*(1 - ne/Ne) + (sigma_IB/ionization_energy.*I - 1/tau_r).*ne;

end