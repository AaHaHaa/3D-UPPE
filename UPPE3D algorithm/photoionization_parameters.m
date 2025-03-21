function s = photoionization_parameters(s)
%PHOTOIONIZATION_PARAMETERS Load material's photoionization parameters

e = 1.60217663e-19; % Coulomb

num_s = length(s.material);
for s_i = 1:num_s
    switch s.material{s_i} % from "NIST Chemistry WebBook"
        case 'H2'
            ionization_energy = 15.42593; % eV
            
            l = 0; % quantum number l
            Z = 1; % effective charge

            me = 9.1093837e-31; % kg
            tau_c = []; % s; averaged collision time of electrons
            tau_r = []; % s; recombination time
            Ne = []; % m^(-3); number density of atoms (neutral electrons) = gas number density, so it's not used
        case 'D2'
            ionization_energy = 15.46658; % eV
            
            l = 0; % quantum number l
            Z = 1; % effective charge

            me = 9.1093837e-31; % kg
            tau_c = []; % s; averaged collision time of electrons
            tau_r = []; % s; recombination time
            Ne = []; % m^(-3); number density of atoms (neutral electrons) = gas number density, so it's not used
        case 'N2'
            ionization_energy = 15.581; % eV
            
            % The values below are taken from https://doi.org/10.1016/S0030-4018(99)00113-3
            % Talebpour et al., "Semi-empirical model for the rate of
            % tunnel ionization of N2 and O2 molecule in an intense
            % Ti:sapphire laser pulse," Opt. Comm. 163, 29-32 (1999)
            l = 0; % quantum number l
            Z = 0.9; % effective charge

            me = 9.1093837e-31; % kg
            tau_c = []; % s; averaged collision time of electrons
            tau_r = []; % s; recombination time
            Ne = []; % m^(-3); number density of atoms (neutral electrons) = gas number density, so it's not used
        case 'O2'
            ionization_energy = 12.0697; % eV
            
            % The values below are taken from https://doi.org/10.1016/S0030-4018(99)00113-3
            % Talebpour et al., "Semi-empirical model for the rate of
            % tunnel ionization of N2 and O2 molecule in an intense
            % Ti:sapphire laser pulse," Opt. Comm. 163, 29-32 (1999)
            l = 0; % quantum number l
            Z = 0.53; % effective charge

            me = 9.1093837e-31; % kg
            tau_c = []; % s; averaged collision time of electrons
            tau_r = []; % s; recombination time
            Ne = []; % m^(-3); number density of atoms (neutral electrons) = gas number density, so it's not used
        case 'CH4'
            ionization_energy = 12.61; % eV
            
            % The ionization of CH4 isn't supported yet.
            l = [];
            Z = [];

            me = 9.1093837e-31; % kg
            tau_c = []; % s; averaged collision time of electrons
            tau_r = []; % s; recombination time
            Ne = []; % m^(-3); number density of atoms (neutral electrons) = gas number density, so it's not used
        case 'He'
            ionization_energy = 24.58741; % eV
            
            l = 0; % quantum number l
            Z = 1; % effective charge

            me = 9.1093837e-31; % kg
            tau_c = []; % s; averaged collision time of electrons
            tau_r = []; % s; recombination time
            Ne = []; % m^(-3); number density of atoms (neutral electrons) = gas number density, so it's not used
        case 'Ne'
            ionization_energy = 21.56454; % eV
            
            l = 1; % quantum number l
            Z = 1; % effective charge

            me = 9.1093837e-31; % kg
            tau_c = []; % s; averaged collision time of electrons
            tau_r = []; % s; recombination time
            Ne = []; % m^(-3); number density of atoms (neutral electrons) = gas number density, so it's not used
        case 'Ar'
            ionization_energy = 15.759; % eV
            
            l = 1; % quantum number l
            Z = 1; % effective charge

            me = 9.1093837e-31; % kg
            tau_c = []; % s; averaged collision time of electrons
            tau_r = []; % s; recombination time
            Ne = []; % m^(-3); number density of atoms (neutral electrons) = gas number density, so it's not used
        case 'Kr'
            ionization_energy = 13.99961; % eV
            
            l = 1; % quantum number l
            Z = 1; % effective charge

            me = 9.1093837e-31; % kg
            tau_c = []; % s; averaged collision time of electrons
            tau_r = []; % s; recombination time
            Ne = []; % m^(-3); number density of atoms (neutral electrons) = gas number density, so it's not used
        case 'Xe'
            ionization_energy = 12.12987; % eV
            
            l = 1; % quantum number l
            Z = 1; % effective charge

            me = 9.1093837e-31; % kg
            tau_c = []; % s; averaged collision time of electrons
            tau_r = []; % s; recombination time
            Ne = []; % m^(-3); number density of atoms (neutral electrons) = gas number density, so it's not used
        case 'CO2'
            ionization_energy = 13.78; % eV
            
            % The ionization of CO2 isn't fully supported yet.
            l = 1; % quantum number l
            Z = 1; % effective charge

            me = 9.1093837e-31; % kg
            tau_c = []; % s; averaged collision time of electrons
            tau_r = []; % s; recombination time
            Ne = []; % m^(-3); number density of atoms (neutral electrons) = gas number density, so it's not used
        case 'N2O'
            ionization_energy = 12.89; % eV
            
            % Assume the values of N2
            % (Copied from N2's)
            l = 0; % quantum number l
            Z = 0.9; % effective charge

            me = 9.1093837e-31; % kg
            tau_c = []; % s; averaged collision time of electrons
            tau_r = []; % s; recombination time
            Ne = []; % m^(-3); number density of atoms (neutral electrons) = gas number density, so it's not used
        case 'silica'
            ionization_energy = 9; % eV
            
            l = []; % quantum number l
            Z = []; % effective charge

            me = 0.86*9.1093837e-31; % kg; silica has a smaller effective electron mass
            tau_c = 1e-15; % s; averaged collision time of electrons
            tau_r = 150e-15; % s; recombination time
            Ne = 2.20e28; % m^(-3); number density of atoms (neutral electrons)
        otherwise
            error('photoionization_parameters:MaterialError',...
                  'This code doesn''t support the ionization computation of the input materials yet');
    end

    s.(s.material{s_i}).ionization = struct('energy',ionization_energy*e,... % J
                                                 'l',l,... % quantum number l
                                                 'Z',Z,... % effective charge
                                                'me',me,... % electron mass
                                                'tau_c',tau_c,... % collision time
                                                'tau_r',tau_r,... % recombination time
                                                'Ne',Ne); % number density of neutral electrons
end

end

