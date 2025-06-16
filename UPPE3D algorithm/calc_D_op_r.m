function [D_op,W_op,loss_op,...
          kc,k0,...
          n,n_pulse,...
          sim] = calc_D_op_r(sim,...
                             n,...
                             Nt,dt,...
                             kr,...
                             Omega,omega,...
                             E_tr,...
                             F_op)
%CALC_D_OP_R It computes the dispersion operator used in 3D-UPPE

E_wr = F_op.Ff(E_tr,[]);

% The dispersion term in the 3D-UPPE, in frequency space
n = ifftshift(n,1); % (Nt,Nr,1,Np); in the order that the fft gives
c = 299792458e-12; % m/ps
k0 = omega/c; % 1/m
k = real(n).*k0; % consider only the refractive index (real part of "n")
% Refractive index is separated into the space-invariant and variant parts.
% The invariant part is taken as the averaged index experienced by the
% field.
if size(n,2) == 1 % n is space-invariant only
    nc = n;
else
    nc = sum(abs(E_wr).^2.*n,2)./sum(abs(E_wr).^2,2);
end
kc = real(nc).*k0;

% Below I set all indices smaller than 1 to 1.
% This is because in solid, refractive index drops around the loss peak,
% down to zero. Both the waveguide operator and Kz_correction (applied
% later) will perform weird, such as having NaN. It is important to neglect
% these indices.
kc(real(nc)<1) = k0(real(nc)<1);

kW2 = k.^2 - kc.^2; % waveguide contribution
W_op = 1i*kW2/2./kc;

% Obtain the omega0 of the input pulse
fftshift_omegas = fftshift(Omega,1);
spectrum = sum(abs(fftshift(E_wr,1)).^2,[2,4]);
omega0 = sum(fftshift_omegas.*spectrum)/sum(spectrum); % 2*pi*THz; the pulse center frequency (under shifted omega)

if ~isfield(sim,'betas')
    if Nt == 1 || ~any(E_wr(:)) % CW case
        sim.betas = [kc;0];
    else
        sim.betas = zeros(2,1,'gpuArray');

        % Obtain the betas of the input pulse
        omega_range = 2*pi/dt; % 2*pi*THz
        omegas_idx_near_pulse = fftshift_omegas>omega0-omega_range/5 & fftshift_omegas<omega0+omega_range/5;% pick only the data near the pulse center frequency to find its beta0 and beta1
        clear spectrum omega0 omega_range;

        fftshift_kc = fftshift(kc,1);
        fit_order = max(2,min(7,sum(omegas_idx_near_pulse)-1)); % 2~7
        [betas_Taylor_coeff,~,mu] = polyfit(fftshift_omegas(omegas_idx_near_pulse),real(fftshift_kc(omegas_idx_near_pulse)),fit_order);
        sim.betas = [betas_Taylor_coeff(end);betas_Taylor_coeff(end-1)];
        sim.betas = [sim.betas(1)-sim.betas(2)*mu(1)/mu(2);...
                     sim.betas(2)/mu(2)];
        clearvars fftshift_omegas fftshift_kc fit_order betas_Taylor_coeff mu
    end
end

Kz = sqrt(complex(kc.^2 - kr.^2)); % the dispersion term; in k-space
D_op = 1i*(Kz-(sim.betas(1)+sim.betas(2)*Omega));

% Remove the high spatial frequencies and fill it with zeros
% It's lossy (negative real part) here, so making it zero is equivalent but
% it makes the adaptive step-size control simpler and faster.
D_op(repmat(kc,1,length(kr),1,1) < kr) = 0;

% In this 3D-UPPE, it is crucial to maintain energy conservation during the
% internal RK4IP for dispersion and waveguide operations. Therefore, the
% operation of loss should be separately computed.
D_op = 1i*imag(D_op);
loss_op = -imag(nc).*k0;

% For photoionization in solids
n_pulse = sum(abs(E_wr).^2.*n,1)./sum(abs(E_wr).^2,1);

end