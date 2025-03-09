function foutput = UPPE3D_propagate_r_gas(fiber, initial_condition, sim, gas)
%UPPE3D_PROPAGATE_R_GAS

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

%% Some Parameters
% Get the numerical parameters from the initial condition.
[Nt,Nr,Np] = size(initial_condition.field,[1,2,4]);
dt = initial_condition.dt;

% During the nonlinear computation, the frequency window is extended to
% three times to prevent aliasing from Raman or FWM. It's shrinked back to
% the original size after the operation.
if Nt == 1 % CW case
    gas_Nt = 1; % = Nt
    gas_dt = dt;
else
    gas_Nt = Nt*3;
    gas_dt = dt/3;
end

%% Check the validity of input parameters
if sim.gpu_yes
    try
        gpuDevice;
    catch
        error('UPPE_propagate_r_gas:GPUError',...
              'No GPU is detected. Please set "sim.gpu_yes=false".');
    end
end

if sim.save_period == 0
    sim.save_period = fiber.L0;
end

num_saves_total = fiber.L0/sim.save_period;
if rem(num_saves_total,1) && rem(num_saves_total+eps(num_saves_total),1) && rem(num_saves_total-eps(num_saves_total),1)
    error('UPPE3D_propagate_r_gas:SizeIncommemsurateError',...
          'The save period is %f m and the fiber length is %f m, which are not commensurate', sim.save_period, fiber.L0)
else
    num_saves_total = round(num_saves_total);
end

if sim.include_Raman
    Fmax = 1/2/dt;
    % Below, only the available Raman (whose preR~=0) are considered.
    % max_Omega and max_T2 are used to determine which Raman model to use (see below).
    switch gas.material
        case 'CH4'
            max_Omega = max(gas.(gas.material).V.omega.*(gas.(gas.material).V.preR~=0));
            max_T2 = max(gas.(gas.material).V.T2.*(gas.(gas.material).V.preR~=0));
        case {'H2','N2','O2'}
            max_Omega = max([gas.(gas.material).R.omega.*(gas.(gas.material).R.preR~=0),...
                             gas.(gas.material).V.omega.*(gas.(gas.material).V.preR~=0)]);
            max_T2 = max([gas.(gas.material).R.T2.*(gas.(gas.material).R.preR~=0),...
                          gas.(gas.material).V.T2.*(gas.(gas.material).V.preR~=0)]);
        case 'air'
            max_Omega = max([gas.N2.R.omega.*(gas.N2.R.preR~=0),...
                             gas.N2.V.omega.*(gas.N2.V.preR~=0),...
                             gas.O2.R.omega.*(gas.O2.R.preR~=0),...
                             gas.O2.V.omega.*(gas.O2.V.preR~=0)]);
            max_T2 = max([gas.N2.R.T2.*(gas.N2.R.preR~=0),...
                          gas.N2.V.T2.*(gas.N2.V.preR~=0),...
                          gas.O2.R.T2.*(gas.O2.R.preR~=0),...
                          gas.O2.V.T2.*(gas.O2.V.preR~=0)]);
        otherwise
            max_Omega = 0;
            max_T2 = 0;
    end
    
    % Check the frequency window
    % It needs to cover 2*f_R to include at least the first-order Stokes 
    % and anti-Stokes Raman scattering.
    if max_Omega/2/pi > Fmax
        error('UPPE_propagate:FrequencyWindowError',...
              ['The frequency window is too small.\n',...
               'It needs to cover 2*f_R to include at least the first-order Stokes and anti-Stokes Raman scattering.']);
    end
    
    % Select the computational model to use
    % 0: Do acyclic convolution for the Raman term; this requires expansion of the time window to prevent aliasing
    % 1: Do cyclic convolution for the Raman term
    %
    % The acyclic-convolution model is implemented because Raman has very 
    % long dephasing time (nanoseconds) where most physics under 
    % investigation happens around fs or ps time scale. Using a large time 
    % window in simulations is a waste of time by 1000 times.
    if Nt*dt < max_T2*10 % 10 is simply some random large number
        gas.model = 0;
    else
        gas.model = 1;
    end
else
    % this is to correctly compute the noise photon in spontaneous_Raman()
    gas.model = 1;
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
               'Fk', @(x,E) FHATHA(x,... % Hankel transform with FHATHA
                                   max(sim.FHATHA.r),...
                                   sim.FHATHA.r,sim.FHATHA.kr,...
                                   sim.FHATHA.dr,sim.FHATHA.dkr,...
                                   sim.FHATHA.l0,sim.FHATHA.exp_prefactor,sim.FHATHA.n2_prefactor,...
                                   sim.FHATHA.ifftQ,...
                                   E),...
              'iFk', @(x,E) FHATHA(x,... % inverse Hankel transform with FHATHA
                                   max(sim.FHATHA.kr),...
                                   sim.FHATHA.kr,sim.FHATHA.r,...
                                   sim.FHATHA.dkr,sim.FHATHA.dr,...
                                   sim.FHATHA.l0,sim.FHATHA.exp_prefactor,sim.FHATHA.n2_prefactor,...
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
 fiber.n,...
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
    fiber.n2 = 0;
    sim.photoionization_model = false;
end
c = 299792458;
permittivity0 = 8.8541878176e-12; % F/m
prefactor_prefactor = (1i/2./kc.*k0.^2/permittivity0);
prefactor = {1 + sim.FHATHA.kr.^2/2./kc.^2,... % correction factor for using kc as the denominator in nonlinear computations (in k-space). Note that all other nonlinear operators are in real space, so it cannot be combined into them here. This is in k-space, so I save it here and will apply it after the nonlinear term is transformed into the k-space.
             prefactor_prefactor.*(permittivity0*2*fiber.n.*ifftshift(fiber.n2,1)/3),... % nonlinear electronic prefactor (in real-xy space)
             prefactor_prefactor./(permittivity0*fiber.n*c)}; % nonlinear Raman prefactor (in real-xy space)

% Photoionization prefactor
if sim.photoionization_model
    me = 9.1093837e-31; % kg; electron mass
    e = 1.60217663e-19; % Coulomb; electron charge 

    prefactor = [prefactor,...
                 {prefactor_prefactor.*(-e^2/me./(omega*1e12).^2),... % ionization-related loss
                  prefactor_prefactor.*(1i*gas.ionization.energy*permittivity0*fiber.n*c./(omega*1e12))}]; % ionization-related phase contribution
end

if sim.gpu_yes
    for i = 1:length(prefactor)
        prefactor{i} = gpuArray(prefactor{i});
    end
end

% Incorporate the damped window to the nonlinear term to remove generation 
% of frequency component near the window edge.
prefactor{2} = prefactor{2}.*sim.damped_window;
prefactor{3} = prefactor{3}.*sim.damped_window;
if sim.photoionization_model
    prefactor{4} = prefactor{4}.*sim.damped_window;
    prefactor{5} = prefactor{5}.*sim.damped_window;
end

%% Zero-padding for upsampling computation of Kerr and Raman effect
% Upsampling to avoid frequency aliasing
%
% Aliasing mostly comes from Raman shift. However, when the spectrum
% becomes broad, aliasing can come from Kerr effect as well due to
% four-wave mixing.
% Instead of applying upsampling to the Raman computation only, it's
% important to apply it to Kerr as well, especially when running
% supercontinuum generation or when your frequency window isn't large
% enough.
[D_op_upsampling,W_op_upsampling,loss_op_upsampling,...
 prefactor,...
 upsampling_zeros] = upsampling(sim,...
                                Nt,gas_Nt,...
                                Nr,Np,...
                                D_op,W_op,loss_op,...
                                prefactor);

%% Set up some parameters for the gas Raman generation equations
if Nt == 1 % ignore Raman scattering under CW cases
    sim.include_Raman = false;
end

% gas_eqn is a container of useful values that will be continuously used during propagation.
[Raw,Rbw,...
 gas_eqn] = Raman_model_for_UPPE_constant_pressure(sim,gas,Nt,...
                                                   gas_Nt,gas_dt,upsampling_zeros,...
                                                   F_op);

%% Setup the exact save points
% We will always save the initial condition as well
save_points = int64(num_saves_total + 1);
save_z = double(0:save_points-1)'*sim.save_period;

%% Photoionization - erfi() lookup table
% Because calculating erfi() is slow, it's faster if I create a lookup table and use interp1().
% The range of input variable for erfi is 0~sqrt(2) only.
if sim.photoionization_model
    n_Am = 10; % the number of summation of Am term in photoionization
    gas_eqn.erfi_x = linspace(0,sqrt(2*(n_Am+1)),1000)';
    gas_eqn.erfi_y = erfi(gas_eqn.erfi_x);
end

%% Modified shot-noise for noise modeling
E_tr_noise_prefactor = shot_noise_r(sim.f0,...
                                    gas_Nt,gas_dt,...
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
 T_delay_out,...
 delta_permittivity,relative_Ne] = SteppingCaller_adaptive_r_gas(sim,gas,gas_eqn,...
                                                                 save_z,save_points,...
                                                                 initial_condition,...
                                                                 prefactor,...
                                                                 F_op,...
                                                                 D_op,...
                                                                 D_op_upsampling, W_op_upsampling, loss_op_upsampling,...
                                                                 Raw, Rbw,...
                                                                 E_tr_noise_prefactor);

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
if sim.include_Raman && sim.scalar
    foutput.delta_permittivity = delta_permittivity;
end
if sim.photoionization_model
    foutput.relative_Ne = relative_Ne;
end

end