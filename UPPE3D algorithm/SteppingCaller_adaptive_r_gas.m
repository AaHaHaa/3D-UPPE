function [E_out,...
          save_z,save_dz,...
          T_delay_out,...
          delta_permittivity,relative_Ne] = SteppingCaller_adaptive_r_gas(sim, gas, gas_eqn,...
                                                                          save_z, save_points,...
                                                                          initial_condition,...
                                                                          prefactor,...
                                                                          F_op,...
                                                                          D_op,...
                                                                          D_op_upsampling, W_op_upsampling, loss_op_upsampling,...
                                                                          Raw, Rbw,...
                                                                          E_tr_noise_prefactor)
%STEPPINGCALLER_ADAPTIVE_R_GAS It starts the pulse propagation.

[Nt,Nr,~,Np] = size(initial_condition.field);
dt = initial_condition.dt;

save_dz = zeros(save_points,1);
T_delay_out = zeros(save_points,1);
if sim.include_Raman && sim.scalar && Nt ~= 1 % non-CW
    delta_permittivity = zeros(Nt,Nr,size(gas_eqn.R_delta_permittivity,2),save_points);
else
    delta_permittivity = []; % dummay variable for output
end
if sim.photoionization_model
    relative_Ne = zeros(Nt,Nr,save_points); % excited electrons due to photoionization
else
    relative_Ne = []; % dummy variable for output
end

% Pulse centering based on the moment of its intensity
if sim.pulse_centering && Nt ~= 1 % non-CW
    % Center the pulse
    temporal_profile = abs(initial_condition.field).^2;
    temporal_profile(temporal_profile<max(temporal_profile,[],1)/10) = 0;
    TCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*temporal_profile,[1,2,4])/sum(temporal_profile(:)));
    % Because circshift is slow on GPU, I discard it.
    if TCenter ~= 0
        if TCenter > 0
            initial_condition.field = [initial_condition.field(1+TCenter:end,:,:,:);initial_condition.field(1:TCenter,:,:,:)];
        elseif TCenter < 0
            initial_condition.field = [initial_condition.field(end+1+TCenter:end,:,:,:);initial_condition.field(1:end+TCenter,:,:,:)];
        end

        if sim.gpu_yes
            TCenter = gather(TCenter);
        end
        T_delay = TCenter*initial_condition.dt;
    else
        T_delay = 0;
    end
else
    T_delay = 0;
end
T_delay_out(1) = T_delay;

E_out = zeros(Nt, Nr, save_points, Np);

% Start by saving the initial condition
if sim.gpu_yes
    E_out(:,:,1,:) = gather(initial_condition.field);
else
    E_out(:,:,1,:) = initial_condition.field;
end

last_E = F_op.Fk(F_op.Ff(initial_condition.field,[]),true); % in k- and frequency space

% Create a progress bar first
if sim.progress_bar
    if ~isfield(sim,'progress_bar_name')
        sim.progress_bar_name = '';
    elseif ~ischar(sim.progress_bar_name)
        error('SteppingCaller_adaptive_r:ProgressBarNameError',...
              '"sim.progress_bar_name" should be a string.');
    end
    h_progress_bar = waitbar(0,sprintf('%s   0.0%%',sim.progress_bar_name),...
        'Name',sprintf('Running 3D-UPPE: %s...',sim.progress_bar_name),...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
    setappdata(h_progress_bar,'canceling',0);

    % Create the cleanup object
    cleanupObj = onCleanup(@()cleanMeUp(h_progress_bar));

    % Use this to control the number of updated time for the progress bar below 1000 times.
    num_progress_updates = 1000;
    progress_bar_z = (1:num_progress_updates)*save_z(end)/num_progress_updates;
    progress_bar_i = 1;
end

z = 0;
save_i = 2; % the 1st one is the initial field
a5 = [];
if ~isfield(sim.adaptive_dz,'max_dz')
    sim.adaptive_dz.max_dz = sim.save_period/10;
end
sim.dz = min(1e-6,sim.adaptive_dz.max_dz); % m; start with a small value to avoid initial blowup
save_dz(1) = sim.dz;
dz_DW = 1e-7; % m; start with a small value to avoid initial blowup
sim.last_dz = 1; % randomly put a number, 1, for initialization

% Then start the propagation
while z+eps(z) < save_z(end) % eps(z) here is necessary due to the numerical error
    % Check for Cancel button press
    if sim.progress_bar && getappdata(h_progress_bar,'canceling')
        error('SteppingCaller_adaptive_r:ProgressBarBreak',...
              'The "cancel" button of the progress bar has been clicked.');
    end

    ever_fail = false;
    previous_E = last_E;
    previous_a5 = a5;

    success = false;
    while ~success
        if ever_fail
            last_E = previous_E;
            a5 = previous_a5;
        end

        [last_E,a5,...
         opt_dz, dz_DW,...
         success] = stepping_RK4IP_adaptive_gas(last_E,...
                                                sim, gas, gas_eqn,...
                                                prefactor,...
                                                F_op,...
                                                D_op_upsampling, W_op_upsampling, loss_op_upsampling,...
                                                Raw, Rbw,...
                                                E_tr_noise_prefactor,...
                                                a5, dz_DW,...
                                                sim.FHATHA.r,...
                                                Nt,dt);

        if ~success
            if opt_dz < 1e-10
                error('SteppingCaller_adaptive_r:adaptiveRK4IPError',...
                      'Adaptive RK4IP continues to fail.\nCheck simulation parameters.');
            end

            ever_fail = true;

            sim.dz = opt_dz;
        end
    end
    sim.last_dz = sim.dz; % previous dz

    % Check for any NaN elements
    if any(isnan(last_E(:))) %any(isnan(last_A),'all')
        error('SteppingCaller_adaptive_r:NaNError',...
              'NaN field encountered, aborting.\nReduce the step size or increase the temporal or frequency window.');
    end

    % Center the pulse
    if sim.pulse_centering && Nt ~= 1 % non-CW
        last_E_in_time = F_op.iFf(last_E,[]);
        temporal_profile = abs(last_E_in_time).^2;
        temporal_profile(temporal_profile<max(temporal_profile,[],1)/10) = 0;
        TCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*temporal_profile,[1,2,4])/sum(temporal_profile(:)));
        % Because circshift is slow on GPU, I discard it.
        if TCenter ~= 0
            if ~isempty(a5) && any(a5(:)) % RK4IP reuses a5 from the previous step
                a5 = F_op.iFf(a5,[]);
            end
            if TCenter > 0
                if ~isempty(a5) && any(a5(:)) % RK4IP reuses a5 from the previous step
                    a5 = F_op.Ff([a5(1+TCenter:end,:,:,:);a5(1:TCenter,:,:,:)],[]);
                end
                last_E = F_op.Ff([last_E_in_time(1+TCenter:end,:,:,:);last_E_in_time(1:TCenter,:,:,:)],[]);
            elseif TCenter < 0
                if ~isempty(a5) && any(a5(:)) % RK4IP reuses a5 from the previous step
                    a5 = F_op.Ff([a5(end+1+TCenter:end,:,:,:);a5(1:end+TCenter,:,:,:)],[]);
                end
                last_E = F_op.Ff([last_E_in_time(end+1+TCenter:end,:,:,:);last_E_in_time(1:end+TCenter,:,:,:)],[]);
            end
            if sim.gpu_yes
                TCenter = gather(TCenter);
            end
            T_delay = T_delay + TCenter*dt;
        end
    end

    % Update z
    z = z + sim.dz;
    % Because the adaptive-step algorithm determines the step size by 
    % checking the error of the spectral intensities from RK3 and RK4, it
    % might ignore changes at the weakest part of the spectrum. This
    % happens in cases of noise-seeded four-wave mixing and noise-seeded
    % Raman scattering far from the pulse central frequency.
    %
    % To account for this effect, I limit the step size to be 10x smaller 
    % than the effective maximum beat length which is
    % 2*pi/(max(eff_betas)-min(eff_betas)).
    % eff_betas is from betas, propagation constants, throughout the time 
    % window but scaled according to the spectral intensity to prevent 
    % taking into account betas where there's no light at all or where 
    % there's some light starting to grow.
    eff_range_D = find_range_D(sum(abs(last_E).^2,2:4),imag(D_op));
    min_beat_length = 2*pi/eff_range_D;
    dz_resolve_beat_length = min_beat_length/4;

    sim.dz = min([opt_dz,save_z(end)-z,sim.adaptive_dz.max_dz,dz_resolve_beat_length]);

    % If it's time to save, get the result from the GPU if necessary,
    % transform to the time domain, and save it
    if z == sim.last_dz
        if sim.include_Raman && sim.scalar && Nt ~= 1 % non-CW
            delta_permittivity(:,:,:,1) = calc_permittivity(sim,gas,gas_eqn,last_E,Nt,F_op);
        end
        if sim.photoionization_model
            E_out_ii = F_op.iFk(F_op.iFf(last_E,[]),true);
            relative_Ne(:,:,1) = calc_Ne(E_out_ii, dt, gas, gas_eqn, sim, F_op);
        end
    end
    if z >= save_z(save_i)-eps(z)
        E_out_ii = F_op.iFk(F_op.iFf(last_E,[]),true);
        if sim.gpu_yes
            save_dz(save_i) = gather(sim.last_dz);
            save_z(save_i) = gather(z);
            E_out(:,:,save_i,:) = gather(E_out_ii);
        else
            save_dz(save_i) = sim.last_dz;
            save_z(save_i) = z;
            E_out(:,:,save_i,:) = E_out_ii;
        end
        if sim.include_Raman && sim.scalar && Nt ~= 1 % non-CW
            delta_permittivity(:,:,:,save_i) = calc_permittivity(sim,gas,gas_eqn,last_E,Nt,F_op);
        end
        if sim.photoionization_model
            relative_Ne(:,:,save_i) = calc_Ne(E_out_ii, dt, gas, gas_eqn, sim, F_op);
        end

        T_delay_out(save_i) = T_delay;

        save_i = save_i + 1;
    end

    % Report current status in the progress bar's message field
    if sim.progress_bar
        if z >= progress_bar_z(progress_bar_i)
            waitbar(gather(z/save_z(end)),h_progress_bar,sprintf('%s%6.1f%%',sim.progress_bar_name,z/save_z(end)*100));
            progress_bar_i = find(z<progress_bar_z,1);
        end
    end
end

end

%% Helper functions
function eff_range_D = find_range_D(spectrum,D)
%FIND_RANGE_D
%
% For an adaptive-dz method, the maximum dz is also limited by the 
% range of the propagation constant, beta0.
% If the FWM, Raman, or anything else happens for multiple frequencies 
% where dz can't resolve their beta0 difference, the outcome can be 
% wrong.
% Here, I multiply the (intensity)^(1/5) of the normalized spectrum to the 
% beta0 to consider beta0 difference of the pulse and exclude those without
% the pulse but within the frequency window. (1/5) is to maximize the 
% contribution of the weak part of the spectrum.

spectrum = spectrum./max(spectrum);

eff_D = D.*spectrum.^(1/5); % I use ^(1/5) to emphasize the weak part
eff_range_D = max(eff_D(:)) - min(eff_D(:));

end
% -------------------------------------------------------------------------
function delta_permittivity = calc_permittivity(sim,gas,gas_eqn, ...
                                                E_wk,...
                                                Nt,...
                                                F_op)
%CALC_PERMITTIVITY finds the Raman-induced permittivity change
%
% Only the imaginary part corresponds to the actual permittiviy contribution of each Raman response.
% The real part is retained so that it's easier to visualize the "intensity" of the phonon strength by taking abs().

E_wr = F_op.iFk(E_wk,true);

num_gas = length(gas.material);
if sim.ellipticity == 0 % linear polarization
    for gas_i = 1:num_gas
        if ismember(gas.material{gas_i},{'H2','N2','O2','air','N2O','CO2'}) % including rotational Raman
            gas_eqn.R_delta_permittivity(:,gas_eqn.cumsum_num_Raman(gas_i)+1) = gas_eqn.R_delta_permittivity(:,gas_eqn.cumsum_num_Raman(gas_i)+1)*4;
        end
    end
end

E_tr_upsampling = F_op.iFf(cat(1,E_wr(1:gas_eqn.n,:,:,:),gas_eqn.upsampling_zeros,E_wr(gas_eqn.n+1:end,:,:,:)),[]);

R_delta_permittivity = permute(gas_eqn.R_delta_permittivity,[1,3,2]); % make it [Nt,1,Raman_type]; The Raman_type dimension are R and V
switch gas.model
    case 0
        delta_permittivity = F_op.iFf(R_delta_permittivity.*ifft(abs(E_tr_upsampling).^2, gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1),[]);
        delta_permittivity = delta_permittivity(gas_eqn.R_downsampling,:,:,:);
    case 1
        delta_permittivity = F_op.iFf(R_delta_permittivity.*F_op.Ff(abs(E_tr_upsampling).^2,[]),[]);
end
% Below follows the computational order as the electric field; however,
% it creates strong aliasing during the second "fft()" operation due to the
% downsampling process for an temporally-elongated delta permittivity.
% To do this correctly, we need to apply the temporal and spectral
% downsampling in the opposite order:
%    1. spectral downsampling
%    2. temporal downsampling with a already-zero-padding delta permittivity
%
% Reverse order creates "spectral downsampling with a non-zero-padding
% delta permittivity" that induces aliasing because delta permittivity
% doesn't reduce to zero at the temporal edges.
% 
% If spectral downsampling is by a factor of an integer, we can simply
% pick data points every "integer" points temporally for spectral
% downsampling. Therefore, only with an integer spectral downsampling
% ratio, both downsampling can be done in the order of temporal-then-
% spectral, downsampling.
%delta_permittivity = ifft(delta_permittivity,[],1).*(permute(max(max(real(sim.mode_profiles.mode_profiles),[],1),[],2),[1,3,2])./mean(sim.mode_profiles.norms,1)).^2; % find the max delta_permittivity of each mode
%delta_permittivity = fft(delta_permittivity([1:gas_eqn.n,gas_eqn.Nt-(Nt-gas_eqn.n-1):gas_eqn.Nt],:,:),[],1); % transform back to time domain

c = 299792458;
permittivity0 = 8.8541878176e-12; % F/m
delta_permittivity = delta_permittivity(1:floor(gas_eqn.Nt/Nt):end,:,:,:)*(2/permittivity0*c); % find the max delta_permittivity

if sim.gpu_yes
    delta_permittivity = gather(delta_permittivity);
end

end
% -------------------------------------------------------------------------
function relative_Ne = calc_Ne(Et, dt, gas, gas_eqn, sim, F_op)

num_gas = length(gas.material);
Ne = 0; % initialization
for gas_i = 1:num_gas
    Ne_i = photoionization_PPT_model(Et, gas.(gas.material{gas_i}).ionization.energy, sim.f0, dt, gas.Ng(gas_i), [],...
                                     gas.(gas.material{gas_i}).ionization.l, gas.(gas.material{gas_i}).ionization.Z,  gas.(gas.material{gas_i}).ionization.me, [],[],...
                                     gas_eqn.erfi_x, gas_eqn.erfi_y,...
                                     sim.ellipticity,...
                                     F_op);
    Ne = Ne + Ne_i;
end
relative_Ne = Ne/sum(gas.Ng);

if sim.gpu_yes
    relative_Ne = gather(relative_Ne);
end

end
% =========================================================================
function cleanMeUp(h_progress_bar)
%CLEANMEUP It deletes the progress bar.

% DELETE the progress bar; don't try to CLOSE it.
delete(h_progress_bar);
    
end