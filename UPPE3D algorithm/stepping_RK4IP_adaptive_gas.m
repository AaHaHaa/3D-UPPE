function [E1,a5,...
          opt_dz,dz_DW,...
          success] = stepping_RK4IP_adaptive_gas(E0,...
                                                 sim, gas, gas_eqn,...
                                                 prefactor,...
                                                 F_op,...
                                                 D_op, W_op, loss_op,...
                                                 Raw, Rbw,...
                                                 E_tr_noise_prefactor,...
                                                 a5_1, dz_DW,...
                                                 r,...
                                                 Nt, dt)
%STEPPING_RK4IP_ADAPTIVE_GAS Take one step with RK4IP without gain

% Noise
E0_tr = F_op.iFk(F_op.iFf(E0,[]),sim.FHATHA.energy_restoration);
I = sum(abs(E0_tr).^2,[1,4]); % size: (1,Nr)
Aeff = (2*pi*trapz(r,I.*r,2))^2/(2*pi*trapz(r,I.^2.*r,2)); % effective mode-field area
E_tr_noise = E_tr_noise_prefactor{1}/sqrt(Aeff);

% Upsampling to avoid frequency aliasing
%
% Aliasing mostly comes from Raman shift. However, when the spectrum
% becomes broad, aliasing can come from Kerr effect as well due to
% four-wave mixing.
% Instead of applying upsampling to the Raman computation only, it's
% important to apply it to Kerr as well, especially when running
% supercontinuum generation or when your frequency window isn't large
% enough.
if Nt ~= 1 % non-CW
    E0_upsampling = cat(1,E0(1:gas_eqn.n,:),gas_eqn.upsampling_zeros,E0(gas_eqn.n+1:end,:));
    if ~isempty(a5_1)
        a5_1_upsampling = cat(1,a5_1(1:gas_eqn.n,:),gas_eqn.upsampling_zeros,a5_1(gas_eqn.n+1:end,:));
    end
else % no upsampling if CW (Nt=1)
    E0_upsampling = E0;
    if ~isempty(a5_1)
        a5_1_upsampling = a5_1;
    end
end

% Represented under the interaction picture
% dz_DW is acquired and used as an initial step for the adaptive-step RK4IP
% for dispersion, diffraction, and waveguide operators
% This speeds up the computation significantly, compared to restarting with
% some random initial guess each time.
[E_IP,...
 dz_DW] = add_DW_RK4IP(E0_upsampling,...
                       F_op,sim.FHATHA.energy_restoration,...
                       D_op,W_op,loss_op,...
                       dz_DW,...
                       sim.dz/2,...
                       prefactor{1},...
                       sim.adaptive_dz.DW_threshold,...
                       r);

% Propagate through the nonlinearity
if isempty(a5_1)
    a5_1_upsampling = N_op(       E0_upsampling,...
                           sim, gas, gas_eqn,...
                           prefactor,...
                           F_op,...
                           Raw, Rbw,...
                           E_tr_noise,...
                           dt/3);
end
a1 = add_DW_RK4IP(a5_1_upsampling,...
                  F_op,sim.FHATHA.energy_restoration,...
                  D_op,W_op,loss_op,...
                  dz_DW,...
                  sim.dz/2,...
                  prefactor{1},...
                  sim.adaptive_dz.DW_threshold,...
                  r);
a2 =       N_op(       E_IP+a1*(sim.dz/2),...
                sim, gas, gas_eqn,...
                prefactor,...
                F_op,...
                Raw, Rbw,...
                E_tr_noise,...
                dt/3);
a3 =       N_op(       E_IP+a2*(sim.dz/2),...
                sim, gas, gas_eqn,...
                prefactor,...
                F_op,...
                Raw, Rbw,...
                E_tr_noise,...
                dt/3);
a4 =       N_op(add_DW_RK4IP(E_IP+a3*(sim.dz),...
                             F_op,sim.FHATHA.energy_restoration,...
                             D_op,W_op,loss_op,...
                             dz_DW,...
                             sim.dz/2,...
                             prefactor{1},...
                             sim.adaptive_dz.DW_threshold,...
                             r),...
                sim, gas, gas_eqn,...
                prefactor,...
                F_op,...
                Raw, Rbw,...
                E_tr_noise,...
                dt/3);

E1 = add_DW_RK4IP(E_IP + (a1+2*a2+2*a3)*(sim.dz/6),...
                  F_op,sim.FHATHA.energy_restoration,...
                  D_op,W_op,loss_op,...
                  dz_DW,...
                  sim.dz/2,...
                  prefactor{1},...
                  sim.adaptive_dz.DW_threshold,...
                  r)...
     + a4*(sim.dz/6);

% Local error estimate
a5 =       N_op(       E1,...
                sim, gas, gas_eqn,...
                prefactor,...
                F_op,...
                Raw, Rbw,...
                E_tr_noise,...
                dt/3);
err = sum(abs((a4(:)-a5(:))*(sim.dz/10)).^2);

% Stepsize control
normE2 = sum(abs(E1(:)).^2);
err = sqrt(err/normE2);
if normE2 == 0 % all-zero field; this will make err NaN (the 2nd if-condition below), so this condition needs to be determined first
    opt_dz = 2*sim.dz;
    success = true;
elseif isnan(normE2) % the computation is just so wrong, so we reduce the step size and do it again
    opt_dz = 0.5*sim.dz;
    success = false;
else
    opt_dz = max(0.5,min(2,0.8*(sim.adaptive_dz.threshold/err)^(1/4)))*sim.dz;

    success = err < sim.adaptive_dz.threshold;
end

% Downsample them back
if Nt ~= 1 % non-CW
    E1 = cat(1,E1(1:gas_eqn.n,:,:,:),E1(end-(Nt-gas_eqn.n-1):end,:,:,:));
    a5 = cat(1,a5(1:gas_eqn.n,:,:,:),a5(end-(Nt-gas_eqn.n-1):end,:,:,:));
end

end

function dEdz_wk = N_op(E_wk,...
                        sim, gas, gas_eqn,...
                        prefactor,...
                        F_op,...
                        Raw, Rbw,...
                        E_tr_noise,...
                        dt)
%N_op Calculate dEdz_wk

E_tr = F_op.iFk(F_op.iFf(E_wk,[]),sim.FHATHA.energy_restoration);
if any(prefactor{2}(:)) % Compute the nonlinearity only when n2 isn't zero
    if size(E_wk,1) ~= 1 % non-CW case
        E_tr_wNoise = E_tr + E_tr_noise;
    
        % Kerr term
        if sim.scalar
            if sim.ellipticity == 0 % linearly polarized
                Kerr =3*E_tr_wNoise.*abs(E_tr_wNoise).^2;
            else % circularly polarized; its Kerr nonlinearity is 2/3 times smaller than the linearly-polarized one
                Kerr =2*E_tr_wNoise.*abs(E_tr_wNoise).^2;
            end
        else
            Kerr = conj(E_tr_wNoise).*sum(E_tr_wNoise.^2,4) + 2*E_tr_wNoise.*sum(abs(E_tr_wNoise).^2,4); % polarization is in the 4th dimension
        end
    
        % Raman term
        if sim.scalar
            Ra = abs(E_tr_wNoise).^2;
        else
            Ra = 2*sum(abs(E_tr_wNoise).^2,4);

            if ~isempty(hbw) % polarized field with an anisotropic Raman
                Rb = E_tr_wNoise.*permute(conj(E_tr_wNoise),[1,2,3,5,4]) + conj(E_tr_wNoise).*permute(E_tr_wNoise,[1,2,3,5,4]);
            end
        end
    
        % Raman term (continued)
        if sim.include_Raman
            switch gas.model
                case 0
                    Ra_upsampling = F_op.Ff(Ra,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt)); % zero-padding for the acyclic convolution theorem to avoid time-domain aliasing
                    Ra_upsampling = cat(1,Ra_upsampling(end-(gas_eqn.m2-1):end,:,:,:,:),Ra_upsampling(1:gas_eqn.m,:,:,:,:)); % only the acyclic_conv_stretch(Nt) points contain useful information; others are from the zero-padding result
                    
                    RaAA = F_op.iFf(Raw.*Ra_upsampling,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt));
                    RaAA = RaAA(gas_eqn.R_downsampling,:,:,:,:).*gas_eqn.phase_shift_acyclic; % select only the valid part
                    
                    if ~sim.scalar % polarized fields with an anisotropic Raman contribution from rotational Raman
                        Rb_upsampling = F_op.Ff(Rb,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt)); % zero-padding for the acyclic convolution theorem to avoid time-domain aliasing
                        Rb_upsampling = cat(1,Rb_upsampling(end-(gas_eqn.m2-1):end,:,:,:,:,:),Rb_upsampling(1:gas_eqn.m,:,:,:,:,:)); % only the acyclic_conv_stretch(Nt) points contain useful information; others are from the zero-padding result
        
                        RbAA = F_op.iFf(Rbw.*Rb_upsampling,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt));
                        RbAA = RbAA(gas_eqn.R_downsampling,:,:,:,:,:).*gas_eqn.phase_shift_acyclic; % select only the valid part
                    end
                case 1
                    Ra = F_op.Ff(Ra,[]); % transform into the frequency domain for upsampling
                    Ra_upsampling = cat(1,Ra(end-(gas_eqn.m2-1):end,:,:,:,:),Ra(1:gas_eqn.n,:,:,:,:));
                    
                    RaAA = F_op.iFf(Raw.*Ra_upsampling,gas_eqn.Nt).*gas_eqn.phase_shift;
                    
                    if ~sim.scalar % polarized fields with an anisotropic Raman contribution from rotational Raman
                        Rb = F_op.Ff(Rb,[]); % transform into the frequency domain for upsampling
                        Rb_upsampling = cat(1,Rb(end-(gas_eqn.m2-1):end,:,:,:,:,:),Rb(1:gas_eqn.n,:,:,:,:,:));
        
                        RbAA = F_op.iFf(Rbw.*Rb_upsampling,gas_eqn.Nt).*gas_eqn.phase_shift;
                    end
            end
            if sim.scalar
                Raman = RaAA.*E_tr_wNoise;
            else % polarized fields with an anisotropic Raman contribution from rotational Raman
                Raman = RaAA.*E_tr_wNoise + sum(RbAA.*permute(E_tr_wNoise,[1,2,3,5,4]),5);
            end
        end
    else % CW case
        % Kerr term
        if sim.scalar
            if sim.ellipticity == 0 % linearly polarized
                Kerr =3*E_tr.*abs(E_tr).^2;
            else % circularly polarized; its Kerr nonlinearity is 2/3 times smaller than the linearly-polarized one
                Kerr =2*E_tr.*abs(E_tr).^2;
            end
        else
            Kerr = conj(E_tr).*sum(E_tr.^2,4) + 2*E_tr.*sum(abs(E_tr).^2,4); % polarization is in the 4th dimension
        end
    
        % Raman term
        if sim.scalar
            Raman = Raw*2*E_tr.*abs(E_tr).^2;
        else
            Raman = Rbw*conj(E_tr).*sum(E_tr.^2,4) + (2*Raw+Rbw)*E_tr.*sum(abs(E_tr).^2,4); % polarization is in the 4th dimension
        end
    end

    if sim.include_Raman
        nonlinear_tr = prefactor{2}.*F_op.Ff(Kerr,[]) + prefactor{3}.*F_op.Ff(Raman,[]);
    else
        nonlinear_tr = prefactor{2}.*F_op.Ff(Kerr,[]);
    end
else
    if sim.gpu_yes
        nonlinear_tr = zeros(size(E_wk),'gpuArray');
    else
        nonlinear_tr = zeros(size(E_wk));
    end
end

if sim.photoionization_model
    inverse_E2 = abs(E_tr).^2;
    inverse_E2(inverse_E2<max(inverse_E2(:))/1e5) = max(inverse_E2(:))/1e5;

    num_gas = length(gas.material);
    nonlinear_photoionization = 0; % initialization
    for gas_i = 1:num_gas
        [Ne,DNeDt] = gas_photoionization_PPT_model(E_tr, gas.(gas.material{gas_i}).ionization.energy, sim.f0, dt, gas.Ng(gas_i),...
                                                   gas.(gas.material{gas_i}).ionization.l, gas.(gas.material{gas_i}).ionization.Z, gas.(gas.material{gas_i}).ionization.me,...
                                                   gas_eqn.erfi_x, gas_eqn.erfi_y,...
                                                   sim.ellipticity,...
                                                   F_op);

        nonlinear_photoionization = nonlinear_photoionization + ...
                                    prefactor{4}.*F_op.Ff(Ne.*E_tr,[]) + gas.(gas.material{gas_i}).ionization.energy*prefactor{5}.*F_op.Ff(DNeDt./inverse_E2.*E_tr,[]);
    end
else
    nonlinear_photoionization = 0;
end

% Finish adding the prefactor
dEdz_wk = prefactor{1}.*F_op.Fk(nonlinear_tr + nonlinear_photoionization,sim.FHATHA.energy_restoration); % nonlinear polarization

end