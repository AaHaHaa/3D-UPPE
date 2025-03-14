function [ne,DneDt] = photoionization_PPT_model(Et, ionization_energy, f0, dt, Ng, n_pulse,...
                                                l, Z, me, tau_c, tau_r,...
                                                erfi_x, erfi_y,...
                                                ellipticity,...
                                                F_op)
%PHOTOIONIZATION_PPT_MODEL This code computes the ionized free electrons
%with the Perelomov-Popov-Terent'ev (PPT) photoionization model. It's used
%as a helper function in the gas-UPPE model.
%   
% Input:
%
%   Et: (Nt,...); the electric field in the time domain
%   ionization_energy: scalar (J)
%   f0: center frequency of the frequency window
%   dt: scalar (ps)
%   Ng: the gas number density (1/m^3)
%   l: quantum number l
%   Z: effective charge
%   erfi_x, erfi_y: (Nt,1); the lookup table for the imaginary error function, erfi()
%   ellipticity: a scalar; only ellipticity=0 is allowed
%
% Output:
%   ne: (Nt,1); the free electron number density (1/m^3)
%   DneDt: (Nt,1); the time derivative of ne (1/m^3/s)

% The polarization dimension is in the 4th in radially-symmetric scheme and
% in 5th dimension in the xy scheme. Here, I kinda try to get them when I
% don't know what exactly the scheme the code is using.
[Nt,Npr,Npxy] = size(Et,[1,4,5]);
if Npr ~= 1 || Npxy ~= 1 || ellipticity ~= 0
    error('Photoionization model works only for linearly polarized fields.');
end
if isempty(l) || ~ismember(l,[0,1]) % l is implemented with 0 or 1 for now
    error('Current material isn''t supported yet.');
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
pulse_phase = real(pulse_phase(floor(Nk/2) : end-floor(Nk/2),:));
% conv() works for only a single-column vector
%pulse_phase = conv(pulse_phase,ones(floor(Nt/100),1)/floor(Nt/100),'same'); 
omega_pulse = -(pulse_phase(3:end,:)-pulse_phase(1:end-2,:))/(2*dt)+2*pi*f0; % THz; I use "central difference" to calculate the slope here
omega_pulse = cat(1,omega_pulse(1,:),omega_pulse,omega_pulse(end,:))*1e12; % Hz

e = 1.60217663e-19; % Coulomb
permittivity0 = 8.85418782e-12; % m^(-3)/kg*s^4*A^2
c = 299792458; % m/s
k = 4*pi*permittivity0;
h = 6.62607015e-34; % m^2*kg/s
hbar = h/2/pi;
a0 = k*hbar^2/me/e^2; % Bohr radius
U_H = e^2/k/a0/2; % hydrogen ionization energy = 13.6 eV

nstar = Z*sqrt(U_H/ionization_energy); % effective principal quantum number n
lstar = max(0,nstar-1); % effective principal quantum number l

I = abs(Et).^2; % intensity; W/m^2

% This modification is for I too small to avoid spurious result in 1./I computation
I(I<max(I(:))/1e5) = max(I(:))/1e5;

ponderomotive_energy = e^2/2/me/permittivity0/c*I./omega_pulse.^2; % because fiber.n is frequency-dependent and is just around one in gases, we ignore it here in the temporal computation
Keldysh_parameter = sqrt(ionization_energy/2./ponderomotive_energy);

kappa = 4*ionization_energy*sqrt(2*me*ionization_energy)/hbar/e;

v = ionization_energy/hbar./omega_pulse.*(1+1/2./Keldysh_parameter.^2);
n_v = ceil(v)-v + permute(0:10,[1,3,2]); % the minimum positive number of n-v+S, where S, an integer, adds n-v until it becomes a positive number

beta = 2*Keldysh_parameter./sqrt(1+Keldysh_parameter.^2);
g = 3/2./Keldysh_parameter.*((1+1/2./Keldysh_parameter.^2).*asinh(Keldysh_parameter) - sqrt(1+Keldysh_parameter.^2)/2./Keldysh_parameter);

% the PPT correction factor
A = cell(1,2); % A = {A0,A1} a cell container for the PPT correction factors, A0 and A1
erfix = interp1(erfi_x,erfi_y,sqrt(beta.*n_v));
A{1} = sum(2/sqrt(3)*Keldysh_parameter.^2./(1+Keldysh_parameter.^2).*exp(-2*n_v.*asinh(Keldysh_parameter)).*erfix,3);
A{1}(Keldysh_parameter<0.8) = 1; % A0 should be close to 1 at small Keldysh parameter
if l == 1
    A{2} = sum(4/sqrt(3*pi)*Keldysh_parameter.^2./(1+Keldysh_parameter.^2).*exp(-2*n_v.*asinh(Keldysh_parameter)).*...
              ( sqrt(pi)/2*beta.*n_v.*erfix + ...
                sqrt(pi)/4*erfix - ...
                sqrt(beta.*n_v)/2.*exp(beta.*n_v) ),3);
    A{2}(Keldysh_parameter<0.8) = 1; % A1 should be close to 1 at small Keldysh parameter
end

Cnl2 = 2^(2*nstar)/nstar/gamma(nstar+lstar+1)/gamma(nstar-lstar);
if isequal(class(Et),'gpuArray')
    W = zeros(size(I),'gpuArray'); % initialize W for the latter summation of the overall ionization rate including m = -l to l
else
    W = zeros(size(I)); % initialize W for the latter summation of the overall ionization rate including m = -l to l
end
for m = -l:l % l is implemented with 0 or 1 for now
    flm = (2*l+1)*gamma(l+abs(m)+1)/2^(abs(m))/factorial(abs(m))/gamma(l-abs(m)+1);
    W = W + Cnl2*flm*ionization_energy/hbar*sqrt(6/pi)*A{abs(m)+1}.*(sqrt(permittivity0*c/2./I)*kappa./sqrt(1+Keldysh_parameter.^2)).^(2*nstar-abs(m)-1.5).*exp(-sqrt(permittivity0*c/2./I)*kappa/3.*g); % because fiber.n is frequency-dependent and is just around one in gases, we ignore it here in the temporal computation
end
W(I<max(I(:))/1e4) = 0; % avoid non-negligible W when there is no or weak field due to numerical precision error

if ~isempty(tau_c) % if solid
    wtc = omega_pulse*tau_c;
    sigma_IB = e^2./n_pulse./omega_pulse/c/permittivity0/me.*wtc./(1+wtc.^2);
    P = cumsum((-W + sigma_IB/ionization_energy.*I - 1/tau_r)*(dt*1e-12),1);
else
    P = cumsum((-W)*(dt*1e-12),1);
end
ne = exp(P).*cumsum(W*Ng.*exp(-P),1)*(dt*1e-12);
if max(ne(:))/Ng > 0.3
    error('photoionization_PPT_model:neError',...
          ['Ion density is too high, reaching %.3f%% of total gas density.\n',...
           'Current implementation supports only perturbative ionization.\n',...
           'Please reduce the peak power or anything that reduces the ionization effect.'],max(ne(:))/Ng*100);
end

if ~isempty(tau_c) % if solid
    DneDt = W.*(Ng - ne) + (sigma_IB/ionization_energy.*I - 1/tau_r).*ne;
else
    DneDt = W.*(Ng - ne);
end

end