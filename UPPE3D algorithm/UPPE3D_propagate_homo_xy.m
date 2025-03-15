function foutput = UPPE3D_propagate_homo_xy(fiber, initial_condition, sim, solid)
%UPPE3D_PROPAGATE_HOMO_XY

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
initial_condition.field = initial_condition.field(:,:,:,end,:); % take only the last field; size: (Nt,Nx,Ny,Nz,Np)

%% Some Parameters
% Get the numerical parameters from the initial condition.
[Nt, Nx, Ny, ~, ~] = size(initial_condition.field);

%% Call solid_info() to load gas parameters
if ~solid.info_called
    [fiber,sim,solid] = solid_info(fiber,sim,solid);
end

%% Check the validity of input parameters
if sim.gpu_yes
    try
        gpuDevice;
    catch
        error('UPPE_propagate_xy:GPUError',...
              'No GPU is detected. Please set "sim.gpu_yes=false".');
    end
end

if sim.save_period == 0
    sim.save_period = fiber.L0;
end

num_saves_total = fiber.L0/sim.save_period;
if rem(num_saves_total,1) && rem(num_saves_total+eps(num_saves_total),1) && rem(num_saves_total-eps(num_saves_total),1)
    error('UPPE3D_propagate_xy:SizeIncommemsurateError',...
          'The save period is %f m and the fiber length is %f m, which are not commensurate', sim.save_period, fiber.L0)
else
    num_saves_total = round(num_saves_total);
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
% Frequency space follow a different Fourier-transform convention from the mathematics
% The convention of k space doesn't matter, so I pick the one that match the mathematics to be consistent with kz
F_op = struct( 'Ff', @(x,n) ifft(x,n,1),...
              'iFf', @(x,n) fft(x,n,1),...
               'Fk', @(x,E) fft(fft(x,[],2),[],3),... % E is used only in the radially-symmetric scheme but just a dummy variable here in the xy scheme
              'iFk', @(x,E) ifft(ifft(x,[],2),[],3)); % E is used only in the radially-symmetric scheme but just a dummy variable here in the xy scheme

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

dx = initial_condition.dx;
if isfield(initial_condition,'dy')
    dy = initial_condition.dy;
else
    dy = dx;
end
if sim.gpu_yes
    dx = gpuArray(dx);
    dy = gpuArray(dy);
end

[D_op,W_op,loss_op,...
 kc,k0,kxy2,...
 fiber.n,n_pulse,...
 sim] = calc_D_op_xy(sim,...
                     fiber.n,...
                     Nt,dt,...
                     Nx,dx,...
                     Ny,dy,...
                     Omega,omega,...
                     initial_condition.field,...
                     F_op);

%% Create a damped frequency window to prevent aliasing
sim.damped_window = create_damped_window_xy(Nt,Nx,Ny);

%% Pre-calculate the factor used in 3D-UPPE
if Nt == 1 % CW case
    sim.photoionization_model = false; % photoionization isn't implemented for CW yet
end
c = 299792458;
permittivity0 = 8.8541878176e-12; % F/m
prefactor_prefactor = (1i/2./kc.*k0.^2/permittivity0);
prefactor = {1 + kxy2/2./kc.^2,... % correction factor for using kc as the denominator in nonlinear computations (in k-space)
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

[fiber,haw,hbw] = solid_Raman_model(fiber,sim,Nt,dt);
if ~sim.include_Raman % no Raman
    fr = 0;
else
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
E_tr_noise_prefactor = shot_noise_xy(sim.f0,...
                                     Nt,dt,...
                                     initial_condition.field,...
                                        dx,dy);
if sim.gpu_yes
    E_tr_noise_prefactor{1} = gpuArray(E_tr_noise_prefactor{1});
end

%% Run the step function over each step
run_start = tic;
% -------------------------------------------------------------------------
[E_out,...
 save_z,save_dz,...
 T_delay_out] = SteppingCaller_adaptive_xy(sim,solid,...
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
    dx = gather(dx);
    dy = gather(dy);
    wait(sim.gpuDevice.Device);
end
fulltime = toc(run_start);

%% Save the results in a struct
foutput = struct('z', save_z,...
                 'dz', save_dz,...
                 'field', E_out,...
                 'dt', initial_condition.dt,...
                 'dx', dx,...
                 'dy', dy,...
                 'betas', sim.betas,...
                 'seconds', fulltime,...
                 't_delay', T_delay_out);

end