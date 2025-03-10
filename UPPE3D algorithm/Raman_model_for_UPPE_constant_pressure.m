function [Raw,Rbw,...
          gas_eqn] = Raman_model_for_UPPE_constant_pressure(sim,gas,Nt,...
                                                            gas_Nt,gas_dt,upsampling_zeros,...
                                                            F_op)
%RAMAN_MODEL_FOR_UPPE It computes the parameters required for the Raman
%computation in gases.
%
%   sim: a structure containing
%       sim.include_Raman
%       sim.gpu_yes
%   gas: a structure containing
%       gas.material
%       gas.model
%       gas.(H2, D2, N2, O2, air, CH4, N2O, CO2)
%   time_window: the size of the time window (ps)
%   Nt: the number of time points in the simulation
%   gas_Nt: the number of time points in the simulation
%   gas_dt: the time interval

num_gas = length(gas.material);
if Nt == 1 % CW case
    num_Raman = zeros(1,num_gas);
    R = []; % initialization
    for gas_i = 1:num_gas
        switch gas.material{gas_i}
            case {'H2','D2','N2','O2'} % rotational + vibrational Raman
                R_i = [sum(gas.(gas.material{gas_i}).R.preR.*gas.(gas.material{gas_i}).R.omega./(gas.(gas.material{gas_i}).R.omega.^2+1./gas.(gas.material{gas_i}).R.T2.^2),2), ...
                       sum(gas.(gas.material{gas_i}).V.preR.*gas.(gas.material{gas_i}).V.omega./(gas.(gas.material{gas_i}).V.omega.^2+1./gas.(gas.material{gas_i}).V.T2.^2),2) ]*1e-12;

                num_Raman_i = 2;
            case 'air'
                R_i = [sum(gas.N2.R.preR.*gas.N2.R.omega./(gas.N2.R.omega.^2+1./gas.N2.R.T2.^2),2) + sum(gas.O2.R.preR.*gas.O2.R.omega./(gas.O2.R.omega.^2+1./gas.O2.R.T2.^2),2), ...
                       sum(gas.N2.V.preR.*gas.N2.V.omega./(gas.N2.V.omega.^2+1./gas.N2.V.T2.^2),2) + sum(gas.O2.V.preR.*gas.O2.V.omega./(gas.O2.V.omega.^2+1./gas.O2.V.T2.^2),2) ]*1e-12;

                num_Raman_i = 2;
            case 'CH4' % only vibrational Raman
                R_i = sum(gas.(gas.material{gas_i}).V.preR.*gas.(gas.material{gas_i}).V.omega./(gas.(gas.material{gas_i}).V.omega.^2+1./gas.(gas.material{gas_i}).V.T2.^2),2)*1e-12;

                num_Raman_i = 1;
            case {'N2O','CO2'} % only rotational Raman
                R_i = sum(gas.(gas.material{gas_i}).R.preR.*gas.(gas.material{gas_i}).R.omega./(gas.(gas.material{gas_i}).R.omega.^2+1./gas.(gas.material{gas_i}).R.T2.^2),2)*1e-12;

                num_Raman_i = 1;
            otherwise
                R_i = [];
                num_Raman_i = 0;
        end

        num_Raman(gas_i) = num_Raman_i;
        R = [R,R_i];
    end
    [Raw,Rbw] = iso_aniso_R(gas.material,num_Raman, ...
                            sim.scalar,sim.ellipticity,...
                            R);
    gas_eqn = [];
else % pulsed case
    time_window = gas_Nt*gas_dt;
    
    acyclic_conv_stretch = @(x) 2*x-1;
    
    n = ceil(Nt/2);
    if sim.include_Raman
        T = (0:gas_Nt-1)'*gas_dt; % ps
        if sim.gpu_yes
            T = gpuArray(T);
        end
        switch gas.model
            case 0
                gas_eqn = struct('Nt', gas_Nt,...
                                 'acyclic_conv_stretch',acyclic_conv_stretch,...
                                 'R_downsampling', [true(gas_Nt,1);false(acyclic_conv_stretch(gas_Nt)-gas_Nt,1)],...
                                 'upsampling_zeros', upsampling_zeros,...
                                 'n',n,'m',ceil(acyclic_conv_stretch(Nt)/2));
                if sim.gpu_yes
                    gas_eqn.R_downsampling = gpuArray(gas_eqn.R_downsampling);
                end
    
                num_Raman = zeros(1,length(gas.material));
                R = []; % initialization
                for gas_i = 1:num_gas
                    switch gas.material{gas_i}
                        case {'H2','D2','N2','O2'} % rotational + vibrational Raman
                            R_i = [sum(gas.(gas.material{gas_i}).R.preR.*exp(-T/gas.(gas.material{gas_i}).R.T2).*exp(1i*gas.(gas.material{gas_i}).R.omega.*T),2),...
                                   sum(gas.(gas.material{gas_i}).V.preR.*exp(-T/gas.(gas.material{gas_i}).V.T2).*exp(1i*gas.(gas.material{gas_i}).V.omega.*T),2)]*acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12); % acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12) is due to the DFT-version convolution theorem
    
                            num_Raman_i = 2;
                        case 'air'
                            R_i = [sum(gas.N2.R.preR.*exp(-T/gas.N2.R.T2).*exp(1i*gas.N2.R.omega.*T),2) + sum(gas.O2.R.preR.*exp(-T/gas.O2.R.T2).*exp(1i*gas.O2.R.omega.*T),2),...
                                   sum(gas.N2.V.preR.*exp(-T/gas.N2.V.T2).*exp(1i*gas.N2.V.omega.*T),2) + sum(gas.O2.V.preR.*exp(-T/gas.O2.V.T2).*exp(1i*gas.O2.V.omega.*T),2)]*acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12); % acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12) is due to the DFT-version convolution theorem
    
                            num_Raman_i = 2;
                        case 'CH4' % only vibrational Raman
                            R_i = sum(gas.(gas.material{gas_i}).V.preR.*exp(-T/gas.(gas.material{gas_i}).V.T2).*exp(1i*gas.(gas.material{gas_i}).V.omega.*T),2)*acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12); % acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12) is due to the DFT-version convolution theorem
    
                            num_Raman_i = 1;
                        case {'N2O','CO2'} % only rotational Raman
                            R_i = sum(gas.(gas.material{gas_i}).R.preR.*exp(-T/gas.(gas.material{gas_i}).R.T2).*exp(1i*gas.(gas.material{gas_i}).R.omega.*T),2)*acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12); % acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12) is due to the DFT-version convolution theorem
    
                            num_Raman_i = 1;
                        otherwise
                            R_i = [];
                            num_Raman_i = 0;
                    end
    
                    num_Raman(gas_i) = num_Raman_i;
                    R = [R,R_i];
                end
                R(isnan(R)) = 0; % in case that some T2=0 such that -T/T2 has a 0/0 term (this happens when the gas pressure is zero)
                % zero-padding in time for the acyclic convolution theorem to avoid time-domain aliasing
                % X = ifft(Y,n,dim) returns the n-point inverse Fourier transform of Y by padding Y with trailing zeros along the dimension "dim" to length n.
                gas_eqn.R_delta_permittivity = F_op.Ff(R,acyclic_conv_stretch(gas_Nt)); % Raman-induced permittivity change
                Rw = F_op.Ff(imag(R),acyclic_conv_stretch(gas_Nt)); % Raman response; only "sin()" part matters in Raman computations
            case 1
                gas_eqn = struct('Nt', gas_Nt,...
                                 'upsampling_zeros', upsampling_zeros,...
                                 'n',n);
    
                num_Raman = zeros(1,length(gas.material));
                R = []; % initialization
                for gas_i = 1:num_gas
                    switch gas.material{gas_i}
                        case {'H2','D2','N2','O2'} % rotational + vibrational Raman
                            R_i = [sum(gas.(gas.material{gas_i}).R.preR.*exp(-T/gas.(gas.material{gas_i}).R.T2).*exp(1i*gas.(gas.material{gas_i}).R.omega.*T),2),...
                                   sum(gas.(gas.material{gas_i}).V.preR.*exp(-T/gas.(gas.material{gas_i}).V.T2).*exp(1i*gas.(gas.material{gas_i}).V.omega.*T),2)]*(time_window*1e-12); % (time_window*1e-12) is due to the DFT-version convolution theorem
    
                            num_Raman_i = 2;
                        case 'air'
                            R_i = [sum(gas.N2.R.preR.*exp(-T/gas.N2.R.T2).*exp(1i*gas.N2.R.omega.*T),2) + sum(gas.O2.R.preR.*exp(-T/gas.O2.R.T2).*exp(1i*gas.O2.R.omega.*T),2),...
                                   sum(gas.N2.V.preR.*exp(-T/gas.N2.V.T2).*exp(1i*gas.N2.V.omega.*T),2) + sum(gas.O2.V.preR.*exp(-T/gas.O2.V.T2).*exp(1i*gas.O2.V.omega.*T),2)]*(time_window*1e-12); % (time_window*1e-12) is due to the DFT-version convolution theorem
    
                            num_Raman_i = 2;
                        case 'CH4' % only vibrational Raman
                            R_i = sum(gas.(gas.material{gas_i}).V.preR.*exp(-T/gas.(gas.material{gas_i}).V.T2).*exp(1i*gas.(gas.material{gas_i}).V.omega.*T),2)*(time_window*1e-12); % (time_window*1e-12) is due to the DFT-version convolution theorem
    
                            num_Raman_i = 1;
                        case {'N2O','CO2'} % only rotational Raman
                            R_i = sum(gas.(gas.material{gas_i}).R.preR.*exp(-T/gas.(gas.material{gas_i}).R.T2).*exp(1i*gas.(gas.material{gas_i}).R.omega.*T),2)*(time_window*1e-12); % (time_window*1e-12) is due to the DFT-version convolution theorem
    
                            num_Raman_i = 1;
                        otherwise
                            R_i = [];
                            num_Raman_i = 0;
                    end
    
                    num_Raman(gas_i) = num_Raman_i;
                    R = [R,R_i];
                end
                R(isnan(R)) = 0; % in case that some T2=0 such that -T/T2 has a 0/0 term (this happens when the gas pressure is zero)
                gas_eqn.R_delta_permittivity = ifft(R,[],1); % Raman-induced permittivity change
                Rw = F_op.Ff(imag(R),[]); % Raman response; only "sin()" part matters in Raman computations
        end
        clear T upsampling_zeros R;
        
        % Raman response (under frequency domain for the convolution operation later)
        [Raw,Rbw,...
         cumsum_num_Raman] = iso_aniso_R(gas.material,num_Raman, ...
                                         sim.scalar,sim.ellipticity,...
                                         Rw);
        gas_eqn.cumsum_num_Raman = cumsum_num_Raman; % for outputing the delta_permittivity in scalar situations
        
        % Consider only the important part during upsampling in time (zero-padding in frequencies)
        % Because of zero-padding spectrally, convolution operation in
        % frequencies includes a lot of multiplications of these zeros. They're
        % removed from the Raman response, gas_eqn.R.
        % Because the time window can be much smaller than the Raman-response
        % decay time, acyclic convolution is implemented here. This results in
        % rearrangement of gas_eqn.R and the corresponding phase shift for a
        % correct result, gas_eqn.phase_shift and gas_eqn.phase_shift_acyclic.
        % Rearrangement of gas_eqn.R is due to the use of MATLAB fft(X,n,2)
        % which zero-pads to n points at the trailing edge; otherwise, we need 
        % to zero-pad ourselves and might be slower than MATLAB internal
        % function, fft. However, an extra phase shift needs to be considered.
        % This might sound confusing. Think carefully about this.
        switch gas.model
            case 0
                m2 = acyclic_conv_stretch(Nt)-gas_eqn.m;
                Raw = [Raw(end-(m2-1):end);Raw(1:gas_eqn.m)];
                if ~isempty(Rbw)
                    Rbw = [Rbw(end-(m2-1):end);Rbw(1:gas_eqn.m)];
                end
                gas_eqn.phase_shift_acyclic = exp(1i*2*pi*m2/acyclic_conv_stretch(gas_Nt)*(0:gas_Nt-1)');
                if sim.gpu_yes
                    gas_eqn.phase_shift_acyclic = gpuArray(gas_eqn.phase_shift_acyclic);
                end
            case 1
                m2 = Nt-n;
                Raw = [Raw(end-(m2-1):end);Raw(1:n)];
                if ~isempty(Rbw)
                    Rbw = [Rbw(end-(m2-1):end);Rbw(1:n)];
                end
                gas_eqn.phase_shift = exp(1i*2*pi*m2/gas_Nt*(0:gas_Nt-1)');
                if sim.gpu_yes
                    gas_eqn.phase_shift = gpuArray(gas_eqn.phase_shift);
                end
        end
        gas_eqn.m2 = m2; % record the index for rearrangement which will be useful later
    else % no Raman
        Raw = []; Rbw = [];
        gas_eqn = struct('Nt',gas_Nt,'upsampling_zeros', upsampling_zeros,'n',n);
        clear upsampling_zeros;
    end
end

end

function [Ra,Rb,...
          cumsum_num_Raman] = iso_aniso_R(material,num_Raman, ...
                                          scalar_yes,ellipticity,...
                                          R)

num_gas = length(material);

R_rot = 0; % initialization
R_vib = 0; % initialization
cumsum_num_Raman = [0,cumsum(num_Raman)];
for gas_i = 1:num_gas
    switch material{gas_i}
        case {'H2','D2','N2','O2','air'} % rotational + vibrational Raman
            R_rot_i = R(:,cumsum_num_Raman(gas_i)+1);
            R_vib_i = R(:,cumsum_num_Raman(gas_i)+2);
        case 'CH4' % only vibrational Raman
            R_rot_i = 0;
            R_vib_i = R(:,cumsum_num_Raman(gas_i)+1);
        case {'N2O','CO2'} % only rotational Raman
            R_rot_i = R(:,cumsum_num_Raman(gas_i)+1);
            R_vib_i = 0;
        otherwise
            R_rot_i = 0;
            R_vib_i = 0;
    end

    R_rot = R_rot + R_rot_i;
    R_vib = R_vib + R_vib_i;
end
clear Rw;
Ra = R_vib - 2*R_rot; % isotropic Raman response
Rb = 6*R_rot; % anisotropic Raman response

% For scalar fields, anisotropic Raman, Rbw, is incorporated into Raw.
if scalar_yes
    if ellipticity == 0 % linear polarization
        Ra = Ra + Rb;
    else % circular polarization: its SRb=SRa/2, so the factor 1/2 is included here
        Ra = Ra + Rb/2;
    end
    Rb = []; % unused dummy variable now
end

end