function [E_out,...
          save_z,save_dz,...
          T_delay_out] = SteppingCaller_adaptive_xy(sim,solid,...
                                                    save_z, save_points,...
                                                    initial_condition,...
                                                    prefactor,...
                                                    F_op,...
                                                    D_op, W_op, loss_op,...
                                                    fr, haw, hbw,...
                                                    E_tr_noise_prefactor,...
                                                    n_pulse)
%STEPPINGCALLER_ADAPTIVE_XY It starts the pulse propagation.

[Nt,Nx,Ny,~,Np] = size(initial_condition.field);
dt = initial_condition.dt;

save_dz = zeros(save_points,1);
T_delay_out = zeros(save_points,1);
if sim.photoionization_model
    relative_Ne = zeros(Nt,Nx,Ny,save_points); % excited electrons due to photoionization
else
    relative_Ne = []; % dummy variable for output
end

% Pulse centering based on the moment of its intensity
if sim.pulse_centering && Nt ~= 1 % non-CW
    % Center the pulse
    temporal_profile = abs(initial_condition.field).^2;
    temporal_profile(temporal_profile<max(temporal_profile,[],1)/10) = 0;
    TCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*temporal_profile,[1,2,3,5])/sum(temporal_profile(:)));
    % Because circshift is slow on GPU, I discard it.
    if TCenter ~= 0
        if TCenter > 0
            initial_condition.field = [initial_condition.field(1+TCenter:end,:,:,:,:);initial_condition.field(1:TCenter,:,:,:,:)];
        elseif TCenter < 0
            initial_condition.field = [initial_condition.field(end+1+TCenter:end,:,:,:,:);initial_condition.field(1:end+TCenter,:,:,:,:)];
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

E_out = zeros(Nt, Nx, Ny, save_points, Np);

% Start by saving the initial condition
if sim.gpu_yes
    E_out(:,:,:,1,:) = gather(initial_condition.field);
else
    E_out(:,:,:,1,:) = initial_condition.field;
end

last_E = F_op.Fk(F_op.Ff(initial_condition.field,[])); % in k- and frequency space

% Create a progress bar first
if sim.progress_bar
    if ~isfield(sim,'progress_bar_name')
        sim.progress_bar_name = '';
    elseif ~ischar(sim.progress_bar_name)
        error('SteppingCaller_adaptive_xy:ProgressBarNameError',...
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
        error('SteppingCaller_adaptive_xy:ProgressBarBreak',...
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
         success] = stepping_RK4IP_adaptive(last_E,...
                                            sim, solid, prefactor,...
                                            F_op,...
                                            D_op, W_op, loss_op,...
                                            fr, haw, hbw,...
                                            E_tr_noise_prefactor,...
                                            a5, dz_DW,...
                                            [],...
                                            dt,n_pulse);

        if ~success
            if opt_dz < 1e-10
                error('SteppingCaller_adaptive_xy:adaptiveRK4IPError',...
                      'Adaptive RK4IP continues to fail.\nCheck simulation parameters.');
            end

            ever_fail = true;

            sim.dz = opt_dz;
        end
    end
    sim.last_dz = sim.dz; % previous dz
    
    % Check for any NaN elements
    if any(any(isnan(last_E))) %any(isnan(last_A),'all')
        error('SteppingCaller_adaptive_xy:NaNError',...
              'NaN field encountered, aborting.\nReduce the step size or increase the temporal or frequency window.');
    end

    % Center the pulse
    if sim.pulse_centering && Nt ~= 1 % non-CW
        last_E_in_time = F_op.iFf(last_E,[]);
        temporal_profile = abs(last_E_in_time).^2;
        temporal_profile(temporal_profile<max(temporal_profile,[],1)/10) = 0;
        TCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*temporal_profile,[1,2,3,5])/sum(temporal_profile(:)));
        % Because circshift is slow on GPU, I discard it.
        if TCenter ~= 0
            if ~isempty(a5) && any(a5(:)) % RK4IP reuses a5 from the previous step
                a5 = F_op.iFf(a5,[]);
            end
            if TCenter > 0
                if ~isempty(a5) && any(a5(:)) % RK4IP reuses a5 from the previous step
                    a5 = F_op.Ff([a5(1+TCenter:end,:,:,:,:);a5(1:TCenter,:,:,:,:)],[]);
                end
                last_E = F_op.Ff([last_E_in_time(1+TCenter:end,:,:,:,:);last_E_in_time(1:TCenter,:,:,:,:)],[]);
            elseif TCenter < 0
                if ~isempty(a5) && any(a5(:)) % RK4IP reuses a5 from the previous step
                    a5 = F_op.Ff([a5(end+1+TCenter:end,:,:,:,:);a5(1:end+TCenter,:,:,:,:)],[]);
                end
                last_E = F_op.Ff([last_E_in_time(end+1+TCenter:end,:,:,:,:);last_E_in_time(1:end+TCenter,:,:,:,:)],[]);
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
    eff_range_D = find_range_D(sum(abs(last_E).^2,2:5),imag(D_op));
    min_beat_length = 2*pi/eff_range_D;
    dz_resolve_beat_length = min_beat_length/4;

    sim.dz = min([opt_dz,save_z(end)-z,sim.adaptive_dz.max_dz,dz_resolve_beat_length]);

    % If it's time to save, get the result from the GPU if necessary,
    % transform to the time domain, and save it
    if z == sim.last_dz
        if sim.photoionization_model
            E_out_ii = F_op.iFk(F_op.iFf(last_E,[]),true);
            relative_Ne(:,:,:,1) = calc_Ne(E_out_ii, dt, solid, sim, F_op, n_pulse);
        end
    end
    if z >= save_z(save_i)-eps(z)
        E_out_ii = F_op.iFk(F_op.iFf(last_E,[]));
        if sim.gpu_yes
            save_dz(save_i) = gather(sim.last_dz);
            save_z(save_i) = gather(z);
            E_out(:,:,:,save_i,:) = gather(E_out_ii);
        else
            save_dz(save_i) = sim.last_dz;
            save_z(save_i) = z;
            E_out(:,:,:,save_i,:) = E_out_ii;
        end
        if sim.photoionization_model
            relative_Ne(:,:,:,save_i) = calc_Ne(E_out_ii, dt, solid, sim, F_op, n_pulse);
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
function relative_Ne = calc_Ne(Et, dt, solid, sim, F_op, n_pulse)

Ne = solid_photoionization_PPT_model(Et, solid.(solid.material{1}).ionization.energy, sim.f0, dt, solid.(solid.material{1}).ionization.Ne, n_pulse,...
                                     solid.(solid.material{1}).ionization.me, solid.(solid.material{1}).ionization.tau_c, solid.(solid.material{1}).ionization.tau_r,...
                                     solid.n_Q, solid.ellip_x, solid.ellipK, solid.ellipE, solid.erfi_x, solid.erfi_y,...
                                     sim.ellipticity,...
                                     F_op);
relative_Ne = Ne/sum(solid.(solid.material{1}).ionization.Ne);

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