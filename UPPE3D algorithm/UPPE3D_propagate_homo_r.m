function foutput = UPPE3D_propagate_homo_r(fiber, initial_condition, sim, solid)
%UPPE3D_PROPAGATE_HOMO_R

%%
if ispc
    sep_char = '\';
else % unix
    sep_char = '/';
end

% Load the folder
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep_char);
current_folder = current_path(1:sep_pos(end-1));

sim.cuda_dir_path = [current_folder 'cuda'];

%% Consider only the last field of the initial condition
initial_condition.field = initial_condition.field(:,:,end,:); % take only the last field; size: (Nt,Nr,Nz,Np)

%% Check the validity of input parameters
if sim.gpu_yes
    try
        gpuDevice;
    catch
        error('UPPE_propagate_r:GPUError',...
              'No GPU is detected. Please set "sim.gpu_yes=false".');
    end
end

if sim.save_period == 0
    sim.save_period = fiber.L0;
end

num_saves_total = fiber.L0/sim.save_period;
if rem(num_saves_total,1) && rem(num_saves_total+eps(num_saves_total),1) && rem(num_saves_total-eps(num_saves_total),1)
    error('UPPE3D_propagate_r:SizeIncommemsurateError',...
          'The save period is %f m and the fiber length is %f m, which are not commensurate', sim.save_period, fiber.L0)
else
    num_saves_total = round(num_saves_total);
end

%% Some Parameters
% Get the numerical parameters from the initial condition.
Nt = size(initial_condition.field,1);

%% Call solid_info() to load solid parameters
if ~solid.info_called
    [fiber,sim,solid] = solid_info(fiber,sim,solid);
end

%% Set up the GPU details
if sim.gpu_yes
    try
        sim.gpuDevice.Device = gpuDevice(sim.gpuDevice.Index); % use the specified GPU device
    catch
        error('Please set the GPU you''re going to use by setting "sim.gpuDevice.Index".');
    end
end

%% Fourier-transform operators
if sim.gpu_yes
    sim.FHATHA.r = gpuArray(sim.FHATHA.r);
    sim.FHATHA.kr = gpuArray(sim.FHATHA.kr);
    sim.FHATHA.dr = gpuArray(sim.FHATHA.dr);
    sim.FHATHA.dkr = gpuArray(sim.FHATHA.dkr);
    sim.FHATHA.l0 = gpuArray(sim.FHATHA.l0);
    sim.FHATHA.exp_prefactor = gpuArray(sim.FHATHA.exp_prefactor);
    sim.FHATHA.r2_prefactor = gpuArray(sim.FHATHA.r2_prefactor);
    sim.FHATHA.ifftQ = gpuArray(sim.FHATHA.ifftQ);
end

% Frequency space follow a different Fourier-transform convention from the mathematics
% The convention of k space doesn't matter, so I pick the one that match the mathematics to be consistent with kz
F_op = struct( 'Ff', @(x,n) ifft(x,n,1),...
              'iFf', @(x,n) fft(x,n,1),...
               'Fk', @(x,E) FHATHA(x,... % Hankel transform with FHATHA
                                   max(sim.FHATHA.r),...
                                   sim.FHATHA.r,sim.FHATHA.kr,...
                                   sim.FHATHA.dr,sim.FHATHA.dkr,...
                                   sim.FHATHA.l0,sim.FHATHA.exp_prefactor,sim.FHATHA.r2_prefactor,...
                                   sim.FHATHA.ifftQ,...
                                   E),...
              'iFk', @(x,E) FHATHA(x,... % inverse Hankel transform with FHATHA
                                   max(sim.FHATHA.kr),...
                                   sim.FHATHA.kr,sim.FHATHA.r,...
                                   sim.FHATHA.dkr,sim.FHATHA.dr,...
                                   sim.FHATHA.l0,sim.FHATHA.exp_prefactor,sim.FHATHA.r2_prefactor,...
                                   sim.FHATHA.ifftQ,...
                                   E));

%% Pre-calculate the dispersion term
% The "Omega" here is offset frequency (omega - omega0), omega: true angular frequency
%                                                        omega0: central angular frequency (=2*pi*f0)
if Nt == 1 % CW case
    dt = 0;
    Omega = 0;
    omega = 2*pi*sim.f0;
else
    if sim.gpu_yes
        dt = gpuArray(initial_condition.dt);
    else
        dt = initial_condition.dt;
    end
    Omega = 2*pi*ifftshift(linspace(-floor(Nt/2), floor((Nt-1)/2), Nt),2)'/(Nt*dt); % in 1/ps, in the order that the ifft gives
    omega = Omega + 2*pi*sim.f0;
end

r = initial_condition.r;
if sim.gpu_yes
    r = gpuArray(r);
end

[D_op,W_op,loss_op,...
 kc,k0,...
 fiber.n,n_pulse,...
 sim] = calc_D_op_r(sim,...
                    fiber.n,...
                    Nt,dt,...
                    sim.FHATHA.kr,...
                    Omega,omega,...
                    initial_condition.field,...
                    F_op);

%% Create a damped frequency window to prevent aliasing
sim.damped_window = create_damped_window_r(Nt,r);

%% Pre-calculate the factor used in 3D-UPPE
if Nt == 1 % CW case
    sim.photoionization_model = false; % photoionization isn't implemented for CW yet
end
c = 299792458;
permittivity0 = 8.8541878176e-12; % F/m
prefactor_prefactor = (1i/2./kc.*k0.^2/permittivity0);
prefactor = {1 + sim.FHATHA.kr.^2/2./kc.^2,... % correction factor for using kc as the denominator in nonlinear computations (in k-space)
             prefactor_prefactor.*(permittivity0*2.*fiber.n*fiber.n2/3)}; % nonlinear prefactor (in real-xy space)

% Photoionization prefactor
if sim.photoionization_model
    e = 1.60217663e-19; % Coulomb; electron charge 
    wtc = (omega*1e12)*solid.(solid.material{1}).ionization.tau_c;

    prefactor = [prefactor,...
                 {prefactor_prefactor.*(-e^2/solid.(solid.material{1}).ionization.me./(omega*1e12).^2.*(-1i+wtc).*wtc./(1+wtc.^2)),... % ionization-related loss
                  prefactor_prefactor.*(1i*permittivity0*fiber.n*c./(omega*1e12))}]; % ionization-related phase contribution
end

if sim.gpu_yes
    for i = 1:length(prefactor)
        prefactor{i} = gpuArray(prefactor{i});
    end
end

% Incorporate the damped window to the n2 term to remove generation of
% frequency component near the window edge.
prefactor{2} = prefactor{2}.*sim.damped_window;
if sim.photoionization_model
    prefactor{3} = prefactor{3}.*sim.damped_window;
    prefactor{4} = prefactor{4}.*sim.damped_window;
end

%% Pre-compute the Raman response in frequency space
if Nt == 1 % ignore Raman scattering under CW cases
    sim.include_Raman = false;
end

if ~sim.include_Raman % no Raman
    fr = 0;
    haw = 0; hbw = 0;
else
    [fiber,haw,hbw] = solid_Raman_model(fiber,sim,Nt,dt);
    fr = fiber.fr;
end

%% Setup the exact save points
% We will always save the initial condition as well
save_points = int64(num_saves_total + 1);
save_z = double(0:save_points-1)'*sim.save_period;

%% Photoionization - erfi() lookup table
% Because calculating erfi() is slow, it's faster if I create a lookup table and use interp1().
% The range of input variable for erfi is 0~sqrt(2) only.
if sim.photoionization_model
    solid.n_Q = 10; % the number of summation terms in photoionization's Q
    solid.ellip_x = linspace(0,1-1e-5,1000)';
    [solid.ellipK,solid.ellipE] = ellipke(solid.ellip_x);
    solid.erfi_x = linspace(0,sqrt(solid.n_Q+2),1000)';
    solid.erfi_y = erfi(solid.erfi_x);
end

%% Modified shot-noise for noise modeling
E_tr_noise_prefactor = shot_noise_r(sim.f0,...
                                    Nt,dt,...
                                    initial_condition.field,...
                                       r);
if sim.gpu_yes
    E_tr_noise_prefactor{1} = gpuArray(E_tr_noise_prefactor{1});
end

%% Run the step function over each step
run_start = tic;
% -------------------------------------------------------------------------
[E_out,...
 save_z,save_dz,...
 T_delay_out] = SteppingCaller_adaptive_r(sim,solid,...
                                          save_z,save_points,...
                                          initial_condition,...
                                          prefactor,...
                                          F_op,...
                                          D_op, W_op, loss_op,...
                                          fr, haw, hbw,...
                                          E_tr_noise_prefactor,...
                                          n_pulse);

% -------------------------------------------------------------------------
% Just to get an accurate timing, wait before recording the time
if sim.gpu_yes
    sim.betas = gather(sim.betas);
    r = gather(r);
    wait(sim.gpuDevice.Device);
end
fulltime = toc(run_start);

%% Save the results in a struct
foutput = struct('z', save_z,...
                 'dz', save_dz,...
                 'field', E_out,...
                 'dt', initial_condition.dt,...
                 'r', r,...
                 'betas', sim.betas,...
                 'seconds', fulltime,...
                 't_delay', T_delay_out);

end