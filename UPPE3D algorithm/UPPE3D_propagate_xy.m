function foutput = UPPE3D_propagate_xy(fiber, initial_condition, sim, dummy_input) %#ok
%UPPE3D_PROPAGATE_XY

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

%% Check the validity of input parameters
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
 fiber.n,...
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
if Nt == 1 % ignore nonlinearity under CW cases
    fiber.n2 = 0;
end
prefactor = {1 + kxy2/2./kc.^2,... % correction factor for using kc as the denominator in nonlinear computations (in k-space)
             (1i/2./kc).*(k0.^2*2.*fiber.n*fiber.n2/3)}; % nonlinear prefactor (in real-xy space)

% Incorporate the damped window to the n2 term to remove generation of
% frequency component near the window edge.
prefactor{2} = prefactor{2}.*sim.damped_window;

%% Pre-compute the Raman response in frequency space
if Nt == 1 % ignore Raman scattering under CW cases
    sim.include_Raman = false;
end

if ~isfield(fiber,'material')
    fiber.material = 'silica';
end
[fiber,haw,hbw] = Raman_model(fiber,sim,Nt,dt);
if ~sim.include_Raman % no Raman
    fr = 0;
else
    fr = fiber.fr;
end

%% Setup the exact save points
% We will always save the initial condition as well
save_points = int64(num_saves_total + 1);
save_z = double(0:save_points-1)'*sim.save_period;

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
 T_delay_out] = SteppingCaller_adaptive_xy(sim,...
                                           save_z,save_points,...
                                           initial_condition,...
                                           prefactor,...
                                           F_op,...
                                           D_op, W_op, loss_op,...
                                           fr, haw, hbw,...
                                           E_tr_noise_prefactor);

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