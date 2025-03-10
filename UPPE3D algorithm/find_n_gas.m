function [n_gas,Raman_absorption] = find_n_gas(material,wavelength,eta)
%FIND_N_GAS It finds the refractive index of gases
%
% Input arguments:
%   material: gas material
%   wavelength: (m)
%   eta: gas density (in amagat)
%
% Output arguments:
%   n_gas: refractive index
%   Raman_absorption: this is, in principle, the imag(n_gas)
%
%   If it's ARHCF, n_gas is real-valued while its imaginary part is in
%   Raman_absorption. Otherwise, Raman_absorption is added to n_gas.

%% Refractive index

% Reference:
% 1. Walter G., et el, "On the Dependence of the Refractive Index of Gases on Temperature" (1903)
% 2. Arthur L. Ruoff and Kouros Ghandehari, "THE REFRACTIVE INDEX OF HYDROGEN AS A FUNCTION OF PRESSURE" (1993)

[a,b] = Sellmeier_coefficients(material); % Sellmeier coefficients
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
switch material
    case {'H2','D2'}
        %n_gas = calc_n_H2(wavelength*1e9,sim.cuda_dir_path,gas.wavelength_order);
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        permittivity_r = n_from_Sellmeier(0.120).^2;
        n_gas_120 = sqrt((permittivity_r - 1)*eta + 1); %  % Sellmeier is valid only above 164nm
        n_gas(wavelength<120e-9) = n_gas_120;
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(material,wavelength,eta); % imag_k_gas
        n_gas = n_gas + 1i*Raman_absorption./(2*pi./wavelength);
    case 'O2'
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        permittivity_r = n_from_Sellmeier(0.4).^2;
        n_gas_400 = sqrt((permittivity_r - 1)*eta + 1); %  % Sellmeier is valid only above 164nm
        n_gas(wavelength<120e-9) = n_gas_400;
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(material,wavelength,eta); % imag_k_gas
        n_gas = n_gas + 1i*Raman_absorption./(2*pi./wavelength);
    case {'air','N2'}
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(material,wavelength,eta); % imag_k_gas
        n_gas = n_gas + 1i*Raman_absorption./(2*pi./wavelength);
    case {'Ar','N2O','CO2'}
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas

        Raman_absorption = 0; % imag_k_gas
    case {'Ne','He'}
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
        
        permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas

        Raman_absorption = 0; % imag_k_gas
    case {'Xe','Kr'}
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
        
        permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        permittivity_r = n_from_Sellmeier(0.113).^2;
        n_gas_113 = sqrt((permittivity_r - 1)*eta + 1); % Sellmeier is valid only above ~113nm
        n_gas(wavelength<113e-9) = n_gas_113;

        Raman_absorption = 0; % imag_k_gas
    case 'CH4'
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
        
        permittivity_r = n_from_Sellmeier(wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        % Avoid the singularity at resonances
        idx_resonance = n_gas < 1;
        n_gas(idx_resonance) = 1;

        Raman_absorption = 0; % imag_k_gas
end

end

