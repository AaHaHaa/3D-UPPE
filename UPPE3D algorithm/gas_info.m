function [fiber,sim,gas] = gas_info(fiber,sim,gas,wavelength)
%GAS_INFO It loads parameters related to gases
%   wavelength: in "m"; from large to small

%% Error check
if length(gas.pressure) ~= length(gas.material)
    error('gas_info_constant_pressure:Pressure_MaterialError',...
          'The length of the gas-pressure (numeric) array should be the same as the length of the gas-material (cell) array');
end
num_gas = length(gas.material);

if any(ismember(gas.material,'air'))
    if any(ismember(gas.material,{'N2','O2'}))
        error('gas_info_constant_pressure:GasMaterialError',...
              ['If gas.material contains ''air'', other materials cannot be ''N2'' or ''O2''.\n',...
               'If necessary, just specify N2 and O2 independently, rather than using ''air''.']);
    end
end

if sim.gpu_yes
    try
        gpuDevice;
    catch
        error('gas_info_constant_pressure:GPUError',...
              'No GPU is detected. Please set "sim.gpu_yes=false".');
    end
end

%% Preset
if ispc
    sep_char = '\';
else % unix
    sep_char = '/';
end
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep_char);
upper_folder = current_path(1:sep_pos(end-1));
addpath(fullfile(upper_folder,'Gas absorption spectra'));

c = 299792458; % m/s
k = 1.38064852e-23; % Boltzmann constant (MKS unit)

pressure0 = 1.01325e5; % Pa
temperature0 = 273.15; % 0 degree Celsius
% eta is calculated with the unit, amagat
% Ideal gas law is used here.
eta = zeros(1,num_gas);
for gas_i = 1:num_gas
    eta(gas_i) = gas.pressure(gas_i)/pressure0*temperature0/gas.temperature;
end

%% Raman parameters

% Number density of the gas
gas.Ng = gas.pressure/k/gas.temperature; % m^(-3)

if any(ismember(gas.material,{'H2','D2','N2','O2','air','CH4','N2O','CO2'}))
    gas = gas_Raman_model(gas,eta); % obtain the Raman parameters according to the gas
else % no Raman
    sim.include_Raman = false;
end

%% Propatation constant (updated with Raman parameters)
fiber.n = 1; % initialization
for gas_i = 1:num_gas
    n_gas_i = find_n_gas(gas.material{gas_i},wavelength,eta(gas_i));
    fiber.n = sqrt(1 + (real(fiber.n).^2-1) + (real(n_gas_i).^2-1)) + 1i*(imag(fiber.n)+imag(n_gas_i));
end

%% Nonlinear coefficient
fiber.n2 = gas_n2(gas.material,wavelength,eta);

%% Ionization potential
if sim.photoionization_model
    gas = photoionization_parameters(gas);
end

%%
% There will be some situations where we want to modify, for example, the 
% Raman response functions, so gas_info() will be called before pulse
% propagation. In this case, we don't want the propagation function to call
% gas_info() again and rewrite the desired parameters set beforehand.
gas.info_called = true; % gas_info() has been called

end