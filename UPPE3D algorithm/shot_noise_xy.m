function E_tr_noise_prefactor = shot_noise_xy(f0,...
                                              Nt,dt,...
                                              field,...
                                                 dx,dy)
%SHOT_NOISE_XY It computes the shot noise included in the governing
% equation of 3D-UPPE.
%
% Because the shot noise approach inherently relies on
% mode-resolved formulation, here we compute only the prefactor. The mode
% area will be included during propagation.

field_size = size(field);

h = 6.62607015e-34; % J*s

time_window = Nt*dt; % ps
f = ifftshift((-floor(Nt/2):floor((Nt-1)/2))'/time_window,1); % THz
real_f = (f+f0)*1e12; % Hz

% I use analytical-signal representation for solving GMMNLSE, so the field
% covers only the positive frequencies.
real_f(real_f<0) = 0; % no noise at negative frequencies

noise_amplitude = sqrt(h*real_f/(time_window*1e-12)); % consider only one photon noise

E_tr_noise_prefactor = {fft(noise_amplitude.*randn(field_size).*exp(1i*2*pi*rand(field_size)),[],1),...
                        dx*dy};

end