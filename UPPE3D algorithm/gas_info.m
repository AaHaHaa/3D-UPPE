function [fiber,sim,gas] = gas_info(fiber,sim,gas,wavelength)
%GAS_INFO It loads parameters related to gases
%   wavelength: in "m"; from small to large

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
eta = gas.pressure/pressure0*temperature0/gas.temperature;

%% Raman parameters

% Number density of the gas
gas.Ng = gas.pressure/k/gas.temperature; % m^(-3)

if ismember(gas.material,{'H2','N2','O2','air','CH4'})
    gas = gas_Raman_model(gas,eta); % Obtain the Raman parameters according to the gas
else % no Raman
    sim.include_Raman = false;
end

%% Propatation constant (updated with Raman parameters)
% Reference:
% 1. Walter G., et el, "On the Dependence of the Refractive Index of Gases on Temperature" (1903)
% 2. Arthur L. Ruoff and Kouros Ghandehari, "THE REFRACTIVE INDEX OF HYDROGEN AS A FUNCTION OF PRESSURE" (1993)

[a,b] = Sellmeier_coefficients(gas.material); % Sellmeier coefficients
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
switch gas.material
    case 'H2'
        %n_gas = calc_n_H2(gas_wavelength*1e9,sim.cuda_dir_path,gas.wavelength_order);
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        permittivity_r = n_from_Sellmeier(0.120).^2;
        n_gas_120 = sqrt((permittivity_r - 1)*eta + 1); %  % Sellmeier is valid only above 164nm
        n_gas(wavelength<120e-9) = n_gas_120;
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(gas.material,wavelength,eta);
        n_gas = n_gas + 1i*Raman_absorption./(2*pi./wavelength);
    case 'O2'
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        permittivity_r = n_from_Sellmeier(0.4).^2;
        n_gas_400 = sqrt((permittivity_r - 1)*eta + 1); %  % Sellmeier is valid only above 400nm
        n_gas(wavelength<120e-9) = n_gas_400;
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(gas.material,wavelength,eta);
        n_gas = n_gas + 1i*Raman_absorption./(2*pi./wavelength);
    case {'air','N2'}
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(gas.material,wavelength,eta);
        n_gas = n_gas + 1i*Raman_absorption./(2*pi./wavelength);
    case {'Ar','Ne','He'}
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
        
        permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
    case {'Xe','Kr'}
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
        
        permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        permittivity_r = n_from_Sellmeier(0.113).^2;
        n_gas_113 = sqrt((permittivity_r - 1)*eta + 1); % Sellmeier is valid only above ~113nm
        n_gas(wavelength<113e-9) = n_gas_113;
    case 'CH4'
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
        
        permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        % Avoid the singularity at resonances
        idx_resonance = n_gas < 1;
        n_gas(idx_resonance) = 1;
end
fiber.n = n_gas;

%% Nonlinear coefficient
% Most values come from
% 1. Shelton, "Nonlinear-optical susceptibilities of gases measured at 1064 and 1319 nm," Phys. Rev. A, 42, 2578-2592 (1990)
% 2. Carsten Bree, Ayhan Demircan, and Gunter Steinmeyer, "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
%
% Note that the values from the following paper are wrong by a factor of 10
% in several of the cited numbers; they copied them wrong.
% Borzsonyi et al., "Measurement of pressure dependent nonlinear refractive index of inert gases," Opt. Express 18, 25847-25854 (2010)
% 
% For Ar, the following two papers seem to require around 7 times weaker n2
% to match their experiments. We have done Ar experiments ourselves which
% show that this reducing factor isn't necessary. We suspect that these old
% works didn't do capillary experiments correctly. It's difficult to align
% the beam into the capillary, which I've been experiencing and is a huge 
% pain. Super straight capillary is required. Super high quality beam is
% also required. Without all these, higher-order modes arise, quickly
% deteriorating the beam spatial quality and the polarization extinction
% ratio. Strong polarization coupling in HE modes in a capillary has
% significant influence to nonlinear propagation.
% 1. Sartania et al.,
%    "Generation of 0.1-TW 5-fs optical pulses at a 1-kHz repetition rate," Opt. Lett. 22(20), 1562-1564 (1997)
% 2. Suda et al.,
%    "Generation of sub-10-fs, 5-mJ-optical pulses using a hollow fiber with a pressure gradient," Appl. Phys. Lett. 86, 111116 (2005)
switch gas.material
    case 'H2' % m^2/(W*atm)
              % This value is taken from Wahlstrand, et al., 
              % "Absolute measurement of the ultrafast nonlinear electronic and rovibrational response in H2 and D2" (2015)
              % Its value, after converted into X3, is close to the paper by Belli et al.,
              % "Vacuum-ultraviolet to infrared supercontinuum in hydrogen-filled photonic crystal fiber" Optica (2015)
              % with X3 = 2.206e-26 m^2/V^2 at standard conditions
        n2 = 0.65e-23;
    case 'N2' % m^2/(W*atm)
              % From Jeffrey M. Brown et al.,
              % "Ab initio calculations of the linear and nonlinear susceptibilities of N2, O2, and air in midinfrared laser pulses"
        P_n2 = 14.63e9; % W
        lambda0_n2 = 0.3334e-6; % m
        n2 = P_n2^(-1)./(lambda0_n2^(-2) - wavelength.^(-2));
        n2(wavelength < 400e-9) = P_n2^(-1)./(lambda0_n2^(-2) - (400e-9).^(-2)); % this equation is valid only above 400 nm
        n2(isinf(n2)) = 0; % avoid the singularity at lambda0_n2
    case 'O2' % m^2/(W*atm)
              % From Jeffrey M. Brown et al.,
              % "Ab initio calculations of the linear and nonlinear susceptibilities of N2, O2, and air in midinfrared laser pulses"
        P_n2 = 14.62e9; % W
        lambda0_n2 = 0.3360e-6; % m
        n2 = P_n2^(-1)./(lambda0_n2^(-2) - wavelength.^(-2));
        n2(wavelength < 400e-9) = P_n2^(-1)./(lambda0_n2^(-2) - (400e-9).^(-2)); % this equation is valid only above 400 nm
        n2(isinf(n2)) = 0; % avoid the singularity at lambda0_n2
    case 'air' % Calculate n2 for N2 and O2
               % Add them up with 79% N2 and 21% O2
        P_n2 = 14.63e9; % W
        lambda0_n2 = 0.3334e-6; % m
        n2_N2 = P_n2^(-1)./(lambda0_n2^(-2) - wavelength.^(-2));
        n2_N2(wavelength < 400e-9) = P_n2^(-1)./(lambda0_n2^(-2) - (400e-9).^(-2)); % this equation is valid only above 400 nm
        P_n2 = 14.62e9; % W
        lambda0_n2 = 0.3360e-6; % m
        n2_O2 = P_n2^(-1)./(lambda0_n2^(-2) - wavelength.^(-2));
        n2_O2(wavelength < 400e-9) = P_n2^(-1)./(lambda0_n2^(-2) - (400e-9).^(-2)); % this equation is valid only above 400 nm
        n2 = 0.79*n2_N2 + 0.21*n2_O2; %clear n2_N2 n2_O2
        n2(isinf(n2)) = 0; % avoid the singularity
    case 'Xe' % m^2/(W*atm)
              % From Shu-Zee Alencious Lo, et al.,
              % "Pulse propagation in hollow-core fiber at high-pressure regime: application to compression of tens of ?J pulses and determination of nonlinear refractive index of xenon at 1.03um" Applied Optics (2018)
        n2 = 50.1e-24;
    case 'Ar' % m^2/(W*atm)
              % From Carsten Bree, Ayhan Demircan, and Gunter Steinmeyer,
              % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
        n2 = 7.96e-24;
    case 'Ne' % m^2/(W*atm)
              % From Carsten Bree, Ayhan Demircan, and Gunter Steinmeyer,
              % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
        n2 = 0.85e-24;
    case 'He' % m^2/(W*atm)
              % From Carsten Bree, Ayhan Demircan, and Gunter Steinmeyer,
              % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
        n2 = 0.34e-24;
    case 'Kr' % m^2/(W*atm)
              % From Carsten Bree, Ayhan Demircan, and Gunter Steinmeyer,
              % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
        n2 = 18.9e-24;
    case 'CH4'
        n2 = 3.118e-23; % m^2/(W*atm)
end
fiber.n2 = n2*eta; % include the gas-density effect

%% Ionization potential
if sim.photoionization_model ~= 0
    switch gas.material % from "NIST Chemistry WebBook"
        case 'H2'
            ionization_energy = 15.42593; % eV
            
            l = 0; % quantum number l
            Z = 1; % effective charge
        case 'N2'
            ionization_energy = 15.581; % eV
            
            % The values below are taken from https://doi.org/10.1016/S0030-4018(99)00113-3
            % Talebpour et al., "Semi-empirical model for the rate of
            % tunnel ionization of N2 and O2 molecule in an intense
            % Ti:sapphire laser pulse," Opt. Comm. 163, 29-32 (1999)
            l = 0; % quantum number l
            Z = 0.9; % effective charge
        case 'O2'
            ionization_energy = 12.0697; % eV
            
            % The values below are taken from https://doi.org/10.1016/S0030-4018(99)00113-3
            % Talebpour et al., "Semi-empirical model for the rate of
            % tunnel ionization of N2 and O2 molecule in an intense
            % Ti:sapphire laser pulse," Opt. Comm. 163, 29-32 (1999)
            l = 0; % quantum number l
            Z = 0.53; % effective charge
        case 'CH4'
            ionization_energy = 12.61; % eV
            
            % The ionization of CH4 isn't supported yet.
            l = [];
            Z = [];
        case 'He'
            ionization_energy = 24.58741; % eV
            
            l = 0; % quantum number l
            Z = 1; % effective charge
        case 'Ne'
            ionization_energy = 21.56454; % eV
            
            l = 1; % quantum number l
            Z = 1; % effective charge
        case 'Ar'
            ionization_energy = 15.759; % eV
            
            l = 1; % quantum number l
            Z = 1; % effective charge
        case 'Kr'
            ionization_energy = 13.99961; % eV
            
            l = 1; % quantum number l
            Z = 1; % effective charge
        case 'Xe'
            ionization_energy = 12.12987; % eV
            
            l = 1; % quantum number l
            Z = 1; % effective charge
        otherwise
            error('This code doesn''t support the ionization computation of other materials yet');
    end
    e = 1.60217663e-19; % Coulomb
    gas.ionization = struct('energy',ionization_energy*e,... % J
                            'l',l,... % quantum number l
                            'Z',Z); % effective charge
end

end