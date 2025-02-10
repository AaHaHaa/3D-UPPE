function gas = gas_Raman_model( gas,eta )
%RAMAN_MODEL It calculates the Raman response of several materials.
%   Input:
%       gas: a structure containing
%           gas.material
%           gas.pressure
%           gas.temperature
%       eta: gas density (amagat)
%
%   Output:
%       gas

c = 299792458; % m/s
h = 6.62607015e-34; % J*s
hbar = h/(2*pi); % J*s
k = 1.38064852e-23; % Boltzmann constant (SI unit)
au_polarizability = 1.64878e-41; % F*m^2; atomic unit

switch gas.material
    case 'H2'
        % T1
        T1 = 1.555e3; % ps
        rot.T1 = T1;
        vib.T1  = T1;

        % T2
        rot.T2 = 1e6/(pi*(6.15/eta+114*eta)); % ps
        vib.T2 = 1e6/(pi*(309/eta*(gas.temperature/298)^0.92+(51.8+0.152*(gas.temperature-298)+4.85e-4*(gas.temperature-298)^2)*eta)); % ps

        % polarizability
        % Below are taken from papers. Don't modify them.
        % If adjustment is needed, change "polarizability_calibration" instead.
        rot.gamma = 2.0239*au_polarizability;
        vib.Dalpha = 3.54e-17;
        vib.Dgamma = 2.57e-17;
        
        % calibration factor
        % To match with the experiments of some papers
        rot.polarizability_calibration = 1.2;
        vib.polarizability_calibration = 1.05; % 1.05 is to fit the vibrational Raman gain from William K. Bischel and Mark J. Dyer's "Wavelength dependence of the absolute Raman gain coefficient for the Q(1) transition in H2"
        % Apply the calibration
        rot.gamma = rot.gamma*rot.polarizability_calibration;
        vib.Dalpha = vib.Dalpha*vib.polarizability_calibration;
        vib.Dgamma = vib.Dgamma*vib.polarizability_calibration;
        
        % nuclear spin statistical constant
        gJ = @(J) mod(J,2)*2 + 1; % 1 if even and 3 if odd
        max_J = 7; % only 7 rotational energy levels below the 1st vibrational energy level
        
        % Energy of each rotational state
        rot.B0 = 59; rot.D0 = 1.6e-2; rot.alpha_e = 3.06; % cm^(-1)
        
        % vibrational Raman shift
        f_vib = 4155; % cm^(-1)
        vib.Omega = 2*pi*f_vib*1e2*c*1e-12; % 2*pi*THz
        
        rot_J = 0:max_J-2;
        vib_J = 0:max_J;
        EJ = @(J) (rot.B0*J.*(J+1) - rot.D0*J.^2.*(J+1).^2)*1e2*h*c; % in the unit "J"
        Z = sum(gJ(vib_J).*(2*vib_J+1).*exp(-EJ(vib_J)/k/gas.temperature)); % partition function considering only the ground vibrational state
        rho = @(J) gJ(J).*exp(-EJ(J)/k/gas.temperature)/Z; % population
        
        % frequency shift of each rotational level
        rot.omega = (EJ(rot_J+2)-EJ(rot_J))/hbar*1e-12; % 2*pi*THz
        vib.omega = vib.Omega - rot.alpha_e*vib_J.*(vib_J+1)*1e2*c*1e-12*2*pi; % 2*pi*THz
        
        % Raman response prefactor
        rot.preR = gas.Ng*1/60/hbar*rot.gamma^2*(rho(rot_J)-rho(rot_J+2)).*(rot_J+2).*(rot_J+1)./(2*rot_J+3);
        vib.preR = gas.Ng*vib.Dalpha^2/4.*(2*vib_J+1).*rho(vib_J)./(vib.omega*1e12);
        
        gas.H2 = struct('R', rot,...
                        'V', vib);
    case {'N2','O2','air'}
        if ismember(gas.material,{'N2','air'})
            % T2
            rot.T2 = 1e6/(pi*3570*eta); % ps
            % Current reference has only T2 = 1e6/(pi*22.5)
            % It doesn't lead to consistent temporal feature as in 
            % "High energy and narrow linewidth N2-filled hollow-core fiber
            % laser at 1.4 μm" by Hong et al. (2023)
            % Therefore, I adjust it myself, following the relation of
            % Raman_linewidth = A/eta+B*eta,
            % whose mininum is at T2_constant with the value 22.5 as in the
            % current reference.
            % Since current reference says that 22.5 should be valid as
            % eta<10, I made it at the minimum point where there is a
            % relatively larger region of points around 22.5.
            T2_constant = 0.02;
            vib.T2 = 1e6/(pi*(11.25*T2_constant/eta +11.25/T2_constant*eta)); %1e6/(pi*22.5); % ps

            % polarizability
            % Below are taken from papers. Don't modify them.
            % If adjustment is needed, change "polarizability_calibration" instead.
            rot.gamma = 4.13*au_polarizability;
            vib.Dalpha = 1.80e-17;
            vib.Dgamma = 2.29e-17;
            
            % calibration factor
            % To match with the experiments of some papers
            rot.polarizability_calibration = 1.17;
            vib.polarizability_calibration = 1.1;
            % Apply the calibration
            rot.gamma = rot.gamma*rot.polarizability_calibration;
            vib.Dalpha = vib.Dalpha*vib.polarizability_calibration;
            vib.Dgamma = vib.Dgamma*vib.polarizability_calibration;
            
            % nuclear spin statistical constant
            gJ = @(J) (2-mod(J,2))*3; % 6 if even and 3 if odd
            max_J = 33; % only 33 rotational energy levels below the 1st vibrational energy level
            
            % Energy of each rotational state7
            rot.B0 = 1.98958; rot.D0 = 5.76e-6; rot.alpha_e = 0.01732; % cm^(-1)
            
            % vibrational Raman shift
            f_vib = 2329.9; % cm^(-1)
            vib.Omega = 2*pi*f_vib*1e2*c*1e-12; % 2*pi*THz
            
            rot_J = 0:max_J-2;
            vib_J = 0:max_J;
            EJ = @(J) (rot.B0*J.*(J+1) - rot.D0*J.^2.*(J+1).^2)*1e2*h*c; % in the unit "J"
            Z = sum(gJ(vib_J).*(2*vib_J+1).*exp(-EJ(vib_J)/k/gas.temperature)); % partition function considering only the ground vibrational state
            rho = @(J) gJ(J).*exp(-EJ(J)/k/gas.temperature)/Z; % population
            
            % frequency shift of each rotational level
            rot.omega = (EJ(rot_J+2)-EJ(rot_J))/hbar*1e-12; % 2*pi*THz
            vib.omega = vib.Omega - rot.alpha_e*vib_J.*(vib_J+1)*1e2*c*1e-12*2*pi; % 2*pi*THz

            % Raman response prefactor
            rot.preR = gas.Ng*1/60/hbar*rot.gamma^2*(rho(rot_J)-rho(rot_J+2)).*(rot_J+2).*(rot_J+1)./(2*rot_J+3);
            vib.preR = gas.Ng*vib.Dalpha^2/4.*(2*vib_J+1).*rho(vib_J)./(vib.omega*1e12);

            if isequal(gas.material,'air')
                N2_ratio_in_air = 0.79;
                rot.preR = rot.preR*N2_ratio_in_air;
                vib.preR = vib.preR*N2_ratio_in_air;
            end
            
            gas.N2 = struct('R', rot,...
                            'V', vib);
        end
        if ismember(gas.material,{'O2','air'})
            % T2
            rot.T2 = 1e6/(pi*1701*eta); % ps
            % T2 should be strongly pressure-dependent; however, there is
            % no existing literature data. I modified it as N2 above.
            T2_constant = 0.02;
            vib.T2 = 1e6/(pi*(27*T2_constant/eta +27/T2_constant*eta)); %1e6/(pi*54); % ps

            % polarizability
            % Below are taken from papers. Don't modify them.
            % If adjustment is needed, change "polarizability_calibration" instead.
            rot.gamma = 6.08*au_polarizability;
            vib.Dalpha = 1.42e-17;
            vib.Dgamma = 2.83e-17;
            
            % calibration factor
            % To match with the experiments of some papers
            rot.polarizability_calibration = 1.17;
            vib.polarizability_calibration = 1.1;
            % Apply the calibration
            rot.gamma = rot.gamma*rot.polarizability_calibration;
            vib.Dalpha = vib.Dalpha*vib.polarizability_calibration;
            vib.Dgamma = vib.Dgamma*vib.polarizability_calibration;
            
            % nuclear spin statistical constant
            gJ = @(J) mod(J,2); % 0 if even and 1 if odd
            max_J = 32; % only 32 rotational energy levels below the 1st vibrational energy level
            
            % Energy of each rotational state
            rot.B0 = 1.43765; rot.D0 = 4.86e-6; rot.alpha_e = 0.01593; % cm^(-1)
            
            % vibrational Raman shift
            f_vib = 1556.3; % cm^(-1)
            vib.Omega = 2*pi*f_vib*1e2*c*1e-12; % 2*pi*THz
            
            rot_J = 0:max_J-2;
            vib_J = 0:max_J;
            EJ = @(J) (rot.B0*J.*(J+1) - rot.D0*J.^2.*(J+1).^2)*1e2*h*c; % in the unit "J"
            Z = sum(gJ(vib_J).*(2*vib_J+1).*exp(-EJ(vib_J)/k/gas.temperature)); % partition function considering only the ground vibrational state
            rho = @(J) gJ(J).*exp(-EJ(J)/k/gas.temperature)/Z; % population
            
            % frequency shift of each rotational level
            rot.omega = (EJ(rot_J+2)-EJ(rot_J))/hbar*1e-12; % 2*pi*THz
            vib.omega = vib.Omega - rot.alpha_e*vib_J.*(vib_J+1)*1e2*c*1e-12*2*pi; % 2*pi*THz

            % Raman response prefactor
            rot.preR = gas.Ng*1/60/hbar*rot.gamma^2*(rho(rot_J)-rho(rot_J+2)).*(rot_J+2).*(rot_J+1)./(2*rot_J+3);
            vib.preR = gas.Ng*vib.Dalpha^2/4.*(2*vib_J+1).*rho(vib_J)./(vib.omega*1e12);
            
            if isequal(gas.material,'air')
                O2_ratio_in_air = 0.21;
                rot.preR = rot.preR*O2_ratio_in_air;
                vib.preR = vib.preR*O2_ratio_in_air;
            end
            
            gas.O2 = struct('R', rot,...
                            'V', vib);
        end
    case 'CH4'
        % Due to the symmetry structure of CH4, it has no rotational Raman
        
        % T2
        vib.T2 = 1e6/(pi*(8220+384*eta)); % ps
        
        % vibrational Raman shift
        f_vib = 2917; % cm^(-1)
        vib.omega = 2*pi*f_vib*1e2*c*1e-12; % 2*pi*THz
        
        % From "Gas Phase Raman Intensities: A Review of "Pre-Laser" Data"
        % by W. F. Murphy, W. Holzer, and H. J. Bernstein.
        % The Raman gain can also be found in 
        % "Measurement of Raman Gain Coefficients of Hydrogen, Deuterium, and Methane"
        % by John J. Ottusch and David A. Rockwel
        % polaribility
        vib.Dalpha = 5.72e-17;
        
        % Raman response prefactor
        vib.preR = gas.Ng*vib.Dalpha^2/4./(vib.omega*1e12);
        
        gas.CH4 = struct('V', vib);
end

end