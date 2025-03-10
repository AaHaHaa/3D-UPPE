function E_tr_noise_prefactor = shot_noise_r(f0,...
                                             Nt,dt,...
                                             field,...
                                                r)
%SHOT_NOISE_R It computes the shot noise included in the governing equation
%of 3D-UPPE.
%
% Because the shot noise approach inherently relies on
% mode-resolved formulation, here we compute only the prefactor. The mode
% area will be included during propagation.

field_size = size(field);
field_size(1) = Nt; % in case that Nt and dt are from gas functions which are bigger due to upsampling to avoid gas-Raman aliasing

if Nt ~= 1
    h = 6.62607015e-34; % J*s
    
    time_window = Nt*dt; % ps
    f = ifftshift((-floor(Nt/2):floor((Nt-1)/2))'/time_window,1); % THz
    real_f = (f+f0)*1e12; % Hz
    
    % I use analytical-signal representation for solving GMMNLSE, so the field
    % covers only the positive frequencies.
    real_f(real_f<0) = 0; % no noise at negative frequencies
    
    noise_amplitude = sqrt(h*real_f/(time_window*1e-12)); % consider only one photon noise
    
    E_tr_noise_prefactor = {fft(noise_amplitude.*randn(field_size).*exp(1i*2*pi*rand(field_size)),[],1)};
else % CW case
    % Noise is based on one noise photon per frequency bin, so it's ignored
    % in CW where there is no such thing as frequency bin.
    if isequal(class(field),'gpuArray')
        E_tr_noise_prefactor = {zeros(field_size,'gpuArray')};
    else
        E_tr_noise_prefactor = {zeros(field_size)};
    end
end

end