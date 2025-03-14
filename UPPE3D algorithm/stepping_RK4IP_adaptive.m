function [E1,a5,...
          opt_dz,dz_DW,...
          success] = stepping_RK4IP_adaptive(E0,...
                                             sim, solid, prefactor,...
                                             F_op,...
                                             D_op, W_op, loss_op,...
                                             fr, haw, hbw,...
                                             E_tr_noise_prefactor,...
                                             a5_1, dz_DW,...
                                             r,...
                                             dt,n_pulse)
%STEPPING_RK4IP_ADAPTIVE Take one step with RK4IP

% Noise
E0_tr = F_op.iFk(F_op.iFf(E0,[]),sim.FHATHA.energy_restoration);
if isempty(r) % use the xy scheme
    I = squeeze(sum(abs(E0_tr).^2,[1,5])); % size: (Nx,Ny)
    Aeff = sum(I*E_tr_noise_prefactor{2},[1,2])^2/sum(I.^2*E_tr_noise_prefactor{2},[1,2]); % effective mode-field area
else % use the radially-symmetric scheme
    I = sum(abs(E0_tr).^2,[1,4]); % size: (1,Nr)
    Aeff = (2*pi*trapz(r,I.*r,2))^2/(2*pi*trapz(r,I.^2.*r,2)); % effective mode-field area
end
E_tr_noise = E_tr_noise_prefactor{1}/sqrt(Aeff);

% Represented under the interaction picture
% dz_DW is acquired and used as an initial step for the adaptive-step RK4IP
% for dispersion, diffraction, and waveguide operators
% This speeds up the computation significantly, compared to restarting with
% some random initial guess each time.
[E_IP,...
 dz_DW] = add_DW_RK4IP(E0,...
                       F_op,sim.FHATHA.energy_restoration,...
                       D_op,W_op,loss_op,...
                       dz_DW,...
                       sim.dz/2,...
                       prefactor{1},...
                       sim.adaptive_dz.DW_threshold,...
                       r);

% Propagate through the nonlinearity
if isempty(a5_1)
    a5_1 = N_op(       E0,...
                sim, solid, prefactor,...
                F_op,...
                fr, haw, hbw,...
                E_tr_noise,...
                r,...
                dt,n_pulse);
end
a1 = add_DW_RK4IP(a5_1,...
                  F_op,sim.FHATHA.energy_restoration,...
                  D_op,W_op,loss_op,...
                  dz_DW,...
                  sim.dz/2,...
                  prefactor{1},...
                  sim.adaptive_dz.DW_threshold,...
                  r);
a2 =       N_op(       E_IP+a1*(sim.dz/2),...
                sim, solid, prefactor,...
                F_op,...
                fr, haw, hbw,...
                E_tr_noise,...
                r,...
                dt,n_pulse);
a3 =       N_op(       E_IP+a2*(sim.dz/2),...
                sim, solid, prefactor,...
                F_op,...
                fr, haw, hbw,...
                E_tr_noise,...
                r,...
                dt,n_pulse);
a4 =       N_op(add_DW_RK4IP(E_IP+a3*(sim.dz),...
                             F_op,sim.FHATHA.energy_restoration,...
                             D_op,W_op,loss_op,...
                             dz_DW,...
                             sim.dz/2,...
                             prefactor{1},...
                             sim.adaptive_dz.DW_threshold,...
                             r),...
                sim, solid, prefactor,...
                F_op,...
                fr, haw, hbw,...
                E_tr_noise,...
                r,...
                dt,n_pulse);

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
                sim, solid, prefactor,...
                F_op,...
                fr, haw, hbw,...
                E_tr_noise,...
                r,...
                dt,n_pulse);
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

end

function dEdz_wk = N_op(E_wk,...
                        sim, solid, prefactor,...
                        F_op,...
                        fr, haw, hbw,...
                        E_tr_noise,...
                        r,...
                        dt,n_pulse)
%N_op Calculate dEdz_wk

E_tr = F_op.iFk(F_op.iFf(E_wk,[]),sim.FHATHA.energy_restoration);
if any(prefactor{2}(:)) % Compute the nonlinearity only when n2 isn't zero
    E_tr_wNoise = E_tr + E_tr_noise;

    % Kerr term
    if sim.scalar
        Kerr = (1-fr)*3*E_tr_wNoise.*abs(E_tr_wNoise).^2;
    else
        if isemtpy(r) % use the xy scheme
            Kerr = (1-fr)*(conj(E_tr_wNoise).*sum(E_tr_wNoise.^2,5) + 2*E_tr_wNoise.*sum(abs(E_tr_wNoise).^2,5)); % polarization is in the 5th dimension
        else % use the radially-symmetric scheme
            Kerr = (1-fr)*(conj(E_tr_wNoise).*sum(E_tr_wNoise.^2,5) + 2*E_tr_wNoise.*sum(abs(E_tr_wNoise).^2,4)); % polarization is in the 4th dimension
        end
    end

    % Raman term
    if sim.include_Raman
        if sim.scalar
            Ra = F_op.iFf(haw.*F_op.Ff(abs(E_tr_wNoise).^2,[]),[]);
            
            nonlinear_tr = Kerr + 3*Ra.*E_tr_wNoise;
        else
            if isemtpy(r) % use the xy scheme
                Ra = F_op.iFf(haw.*F_op.Ff(sum(abs(E_tr_wNoise).^2,5),[]),[]);
            else % use the radially-symmetric scheme
                Ra = F_op.iFf(haw.*F_op.Ff(sum(abs(E_tr_wNoise).^2,4),[]),[]);
            end

            if isempty(hbw)
                nonlinear_tr = Kerr + 3*Ra.*E_tr_wNoise;
            else % polarized field with an anisotropic Raman
                if isemtpy(r) % use the xy scheme
                    Rb = F_op.iFf(hbw.*F_op.Ff(E_tr_wNoise.*permute(conj(E_tr_wNoise),[1,2,3,4,6,5]) + conj(E_tr_wNoise).*permute(E_tr_wNoise,[1,2,3,4,6,5]),[]),[]);
    
                    nonlinear_tr = Kerr + ( 3*Ra.*E_tr_wNoise + 3/2*sum(Rb.*permute(E_tr_wNoise,[1,2,3,4,6,5]),6) );
                else % use the radially-symmetric scheme
                    Rb = F_op.iFf(hbw.*F_op.Ff(E_tr_wNoise.*permute(conj(E_tr_wNoise),[1,2,3,5,4]) + conj(E_tr_wNoise).*permute(E_tr_wNoise,[1,2,3,5,4]),[]),[]);
    
                    nonlinear_tr = Kerr + ( 3*Ra.*E_tr_wNoise + 3/2*sum(Rb.*permute(E_tr_wNoise,[1,2,3,5,4]),5) );
                end
            end
        end
    else
        nonlinear_tr = Kerr;
    end
    nonlinear_tr = prefactor{2}.*F_op.Ff(nonlinear_tr,[]);
else
    nonlinear_tr = 0;
end

if sim.photoionization_model
    inverse_E2 = abs(E_tr).^2;
    inverse_E2(inverse_E2<max(inverse_E2(:))/1e5) = max(inverse_E2(:))/1e5;

    [Ne,DNeDt] = solid_photoionization_PPT_model(E_tr, solid.(solid.material{1}).ionization.energy, sim.f0, dt, solid.(solid.material{1}).ionization.Ne, n_pulse,...
                                                 solid.(solid.material{1}).ionization.me, solid.(solid.material{1}).ionization.tau_c, solid.(solid.material{1}).ionization.tau_r,...
                                                 solid.n_Q, solid.ellip_x, solid.ellipK, solid.ellipE, solid.erfi_x, solid.erfi_y,...
                                                 sim.ellipticity,...
                                                 F_op);

    nonlinear_photoionization = prefactor{3}.*F_op.Ff(Ne.*E_tr,[]) + solid.(solid.material{1}).ionization.energy*prefactor{4}.*F_op.Ff(DNeDt./inverse_E2.*E_tr,[]);
else
    nonlinear_photoionization = 0;
end

% Finish adding the prefactor
dEdz_wk = prefactor{1}.*F_op.Fk(nonlinear_tr + nonlinear_photoionization,sim.FHATHA.energy_restoration); % nonlinear polarization

end