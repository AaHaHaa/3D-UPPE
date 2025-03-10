function gas = gas_Raman_model( gas,eta )
%GAS_RAMAN_MODEL It calculates the Raman response of several materials.
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

num_gas = length(gas.material);
for gas_i = 1:num_gas
    switch gas.material{gas_i}
        case 'H2'
            % T1
            T1 = 1.555e3; % ps
            rot.T1 = T1;
            vib.T1  = T1;
    
            % T2
            % From Herring et al., "Temperature and density dependence of the linewidths and line shifts of the rotational Raman lines in N_2 and H_2," Phys. Rev. A 34, 1944-1951 (1986)
            rot.T2 = 1e6/(pi*(6.15/eta(gas_i)+114*eta(gas_i))); % ps
            % From Bischel and Dyer, "Temperature dependence of the Raman linewidth and line shift for the Q(1) and Q(0) transitions in normal and para-H_2," Phys. Rev. A 33, 3113-3123 (1986)
            vib.T2 = 1e6/(pi*(309/eta(gas_i)*(gas.temperature/298)^0.92+(51.8+0.152*(gas.temperature-298)+4.85e-4*(gas.temperature-298)^2)*eta(gas_i))); % ps
    
            % polarizability
            % Below are taken from papers. Don't modify them.
            % If adjustment is needed, change "polarizability_calibration" instead.
            %
            % From Kolos and Wolniewicz, "Polarizability of the Hydrogen Molecule," J. Chem. Phys. 46, 1426-1432 (1967)
            % gamma are directly taken from gamma_0000 in Table III.
            % Dalpha and Dgamma are computed from finding the slope with data in Table II.
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
            % From Hanson and Poirier, "Stimulated rotational Raman conversion in H_2, D_2, and HD," IEEE J. Quantum Electron. 29, 2342-345 (1993)
            gJ = @(J) mod(J,2)*2 + 1; % 1 if even and 3 if odd
            max_J = 7; % only 7 rotational energy levels below the 1st vibrational energy level
            
            % Energy of each rotational state
            % B0 and D0 are from 
            %    Wahlstrand et al., "Absolute measurement of the ultrafast nonlinear electronic and rovibrational response in H_2 and D_2," Phys. Rev. A 92, 063828 (2015)
            % alpha_e is from
            %    W. Demtröder, “Diatomic Molecules,” in Atoms, Molecules and Photons: An Introduction to Atomic, Molecular- and Quantum Physics, (Springer Berlin Heidelberg, Berlin, Heidelberg, 2010), pp. 327-381.
            rot.B0 = 58.9; rot.D0 = 0.05; rot.alpha_e = 3.06; % cm^(-1)
            
            % vibrational Raman shift
            % From Minck et al., "LASER-STIMULATED RAMAN EFFECT AND RESONANT FOUR-PHOTON INTERACTIONS IN GASES H_2, D_2, AND CH_4," Appl. Phys. Lett. 3, 181-184 (1963)
            f_vib = 4155.21; % cm^(-1)
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
            rot.preR = gas.Ng(gas_i)*1/60/hbar*rot.gamma^2*(rho(rot_J)-rho(rot_J+2)).*(rot_J+2).*(rot_J+1)./(2*rot_J+3);
            vib.preR = gas.Ng(gas_i)*vib.Dalpha^2/4.*(2*vib_J+1).*rho(vib_J)./(vib.omega*1e12);
            
            gas.H2 = struct('R', rot,...
                            'V', vib);
        case 'D2'
            % Some note:
            % By plotting rot.omega vs. rot.preR, we can see that S(2)
            % dominates in D2, in contrast to H2's S(1).
            %
            %
            % T2
            % From Hanson and Poirier, "Stimulated rotational Raman conversion in H_2, D_2, and HD," IEEE J. Quantum Electron. 29, 2342-345 (1993)
            % It contains only the 134 factor.
            % From H2's vib.T2 relation, it seems that 51.8 --> 120 is accompanied by 309 --> 101 with approximately the same factor, possibly due to the inverse eta relation.
            % From this speculation, I scale D2's rot.T2 accordingly.
            % H2's 114 --> D2's 134 leads to H2's 6.15 --> D2's 6.15*(134/114)=7.2, so I put 7 below.
            rot.T2 = 1e6/(pi*(7/eta(gas_i)+134*eta(gas_i))); % ps
            vib.T2 = 1e6/(pi*(101/eta(gas_i)+120*eta(gas_i))); % ps
    
            % polarizability
            % Below are taken from papers. Don't modify them.
            % If adjustment is needed, change "polarizability_calibration" instead.
            %
            % From Kolos and Wolniewicz, "Polarizability of the Hydrogen Molecule," J. Chem. Phys. 46, 1426-1432 (1967)
            % gamma are directly taken from gamma_0000 in Table III.
            rot.gamma = 1.9584*au_polarizability;
            % From Wahlstrand et al., "Absolute measurement of the ultrafast nonlinear electronic and rovibrational response in H_2 and D_2," Phys. Rev. A 92, 063828 (2015)
            % Wahlstand has the value of 1.22 for D2 and 1.30 for H2, so I just do linear scaling. Lazy to convert between units.
            vib.Dalpha = 3.54e-17*1.22/1.30;
            % From Raj et al., "Polarizability tensor invariants of H_2, HD, and D_2," J. Chem. Phys. 148, 104308 (2018)
            % Because their units are different what I used here, some conversion based on H2's value I have here needs to be done.
            vib.Dgamma = (2.0254/2.0947)*(2.57e-17/(2.0239*au_polarizability))*1.9584*au_polarizability;
            
            % calibration factor
            % (Copied from H2's.)
            % To match with the experiments of some papers
            rot.polarizability_calibration = 1.2;
            vib.polarizability_calibration = 1.05; % 1.05 is to fit the vibrational Raman gain from William K. Bischel and Mark J. Dyer's "Wavelength dependence of the absolute Raman gain coefficient for the Q(1) transition in H2"
            % Apply the calibration
            rot.gamma = rot.gamma*rot.polarizability_calibration;
            vib.Dalpha = vib.Dalpha*vib.polarizability_calibration;
            vib.Dgamma = vib.Dgamma*vib.polarizability_calibration;
            
            % nuclear spin statistical constant
            % From Hanson and Poirier, "Stimulated rotational Raman conversion in H_2, D_2, and HD," IEEE J. Quantum Electron. 29, 2342-345 (1993)
            gJ = @(J) -mod(J,2)*3 + 6; % 6 if even and 3 if odd
            max_J = 10; % 10 rotational energy levels are good enough; above them, the Raman strength (preR below) is too small
            
            % Energy of each rotational state
            % B0 and D0 are from 
            %    Wahlstrand et al., "Absolute measurement of the ultrafast nonlinear electronic and rovibrational response in H_2 and D_2," Phys. Rev. A 92, 063828 (2015)
            % alpha_e is from
            %    W. Demtröder, “Diatomic Molecules,” in Atoms, Molecules and Photons: An Introduction to Atomic, Molecular- and Quantum Physics, (Springer Berlin Heidelberg, Berlin, Heidelberg, 2010), pp. 327-381.
            %
            % I didnt' find the value for D2's alpha_e, so I set it zero.
            rot.B0 = 29.3; rot.D0 = 0.021; rot.alpha_e = 0; % cm^(-1)
            
            % vibrational Raman shift
            % From Minck et al., "LASER-STIMULATED RAMAN EFFECT AND RESONANT FOUR-PHOTON INTERACTIONS IN GASES H_2, D_2, AND CH_4," Appl. Phys. Lett. 3, 181-184 (1963)
            f_vib = 2991.39; % cm^(-1)
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
            rot.preR = gas.Ng(gas_i)*1/60/hbar*rot.gamma^2*(rho(rot_J)-rho(rot_J+2)).*(rot_J+2).*(rot_J+1)./(2*rot_J+3);
            vib.preR = gas.Ng(gas_i)*vib.Dalpha^2/4.*(2*vib_J+1).*rho(vib_J)./(vib.omega*1e12);
            
            gas.D2 = struct('R', rot,...
                            'V', vib);
        case {'N2','O2','air'}
            if ismember(gas.material{gas_i},{'N2','air'})
                % T2
                rot.T2 = 1e6/(pi*3570*eta(gas_i)); % ps
                % Current reference has only T2 = 1e6/(pi*22.5)
                % It doesn't lead to consistent temporal feature as in 
                % "High energy and narrow linewidth N2-filled hollow-core fiber
                % laser at 1.4 μm" by Hong et al. (2023)
                % Therefore, I adjust it myself, following the relation of
                % Raman_linewidth = A/eta(gas_i)+B*eta(gas_i),
                % whose mininum is at T2_constant with the value 22.5 as in the
                % current reference.
                % Since current reference says that 22.5 should be valid as
                % eta(gas_i)<10, I made it at the minimum point where there is a
                % relatively larger region of points around 22.5.
                T2_constant = 0.02;
                vib.T2 = 1e6/(pi*(11.25*T2_constant/eta(gas_i) +11.25/T2_constant*eta(gas_i))); %1e6/(pi*22.5); % ps
    
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
                rot.preR = gas.Ng(gas_i)*1/60/hbar*rot.gamma^2*(rho(rot_J)-rho(rot_J+2)).*(rot_J+2).*(rot_J+1)./(2*rot_J+3);
                vib.preR = gas.Ng(gas_i)*vib.Dalpha^2/4.*(2*vib_J+1).*rho(vib_J)./(vib.omega*1e12);
    
                if isequal(gas.material{gas_i},'air')
                    N2_ratio_in_air = 0.79;
                    rot.preR = rot.preR*N2_ratio_in_air;
                    vib.preR = vib.preR*N2_ratio_in_air;
                end
                
                gas.N2 = struct('R', rot,...
                                'V', vib);
            end
            if ismember(gas.material{gas_i},{'O2','air'})
                % T2
                rot.T2 = 1e6/(pi*1701*eta(gas_i)); % ps
                % T2 should be strongly pressure-dependent; however, there is
                % no existing literature data. I modified it as N2 above.
                T2_constant = 0.02;
                vib.T2 = 1e6/(pi*(27*T2_constant/eta(gas_i) +27/T2_constant*eta(gas_i))); %1e6/(pi*54); % ps
    
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
                rot.preR = gas.Ng(gas_i)*1/60/hbar*rot.gamma^2*(rho(rot_J)-rho(rot_J+2)).*(rot_J+2).*(rot_J+1)./(2*rot_J+3);
                vib.preR = gas.Ng(gas_i)*vib.Dalpha^2/4.*(2*vib_J+1).*rho(vib_J)./(vib.omega*1e12);
                
                if isequal(gas.material{gas_i},'air')
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
            vib.T2 = 1e6/(pi*(8220+384*eta(gas_i))); % ps
            
            % vibrational Raman shift
            % From Minck et al., "LASER-STIMULATED RAMAN EFFECT AND RESONANT FOUR-PHOTON INTERACTIONS IN GASES H_2, D_2, AND CH_4," Appl. Phys. Lett. 3, 181-184 (1963)
            f_vib = 2916.5; % cm^(-1)
            vib.omega = 2*pi*f_vib*1e2*c*1e-12; % 2*pi*THz
            
            % From "Gas Phase Raman Intensities: A Review of "Pre-Laser" Data"
            % by W. F. Murphy, W. Holzer, and H. J. Bernstein.
            % The Raman gain can also be found in 
            % "Measurement of Raman Gain Coefficients of Hydrogen, Deuterium, and Methane"
            % by John J. Ottusch and David A. Rockwel
            % polaribility
            vib.Dalpha = 5.72e-17;
            
            % Raman response prefactor
            vib.preR = gas.Ng(gas_i)*vib.Dalpha^2/4./(vib.omega*1e12);
            
            gas.CH4 = struct('V', vib);
        case 'N2O'
            % There is a lack of data for N2O, so the dephasing time T2 is
            % copied from N2. Also, vibrational Raman is ignored since I cannot
            % find its data.
            % T2
            rot.T2 = 1e6/(pi*3570*eta(gas_i)); % ps
    
            % polarizability
            % From
            % 1. Dempsey et al., "Comparative study of optical nonlinearitie of CO2 and N2O via single-shot spatio-temporally-resolved visualization," Opt. Commun. 545, 129669 (2023)
            % 2. Truong et al., "Spectral Broadening and Pulse Compression in Molecular Gas-Filled Hollow-Core Fibers," IEEE J. Sel. Topics Quantum Electron. 30(6), 1-11 (2024)
            % I feel lazy to convert values between units, so I obtain N2O's gamma by scaling it with N2's value in the second paper
            N2_paper = 9.3;
            N2O_paper = 22.7;
            rot.gamma = 4.13*au_polarizability*(N2O_paper/N2_paper); % 4.13 is from N2's
            
            % calibration factor
            % To match with the experiments of some papers
            % (I haven't tested yet, so I put one)
            rot.polarizability_calibration = 1;
            % Apply the calibration
            rot.gamma = rot.gamma*rot.polarizability_calibration;
            
            % nuclear spin statistical constant
            gJ = @(J) 1; % always 1
            max_J = 100; % only 33 rotational energy levels below the 1st vibrational energy level
            
            % Energy of each rotational state7
            rot.B0 = 0.41; % cm^(-1)
            
            rot_J = 0:max_J-2;
            vib_J = 0:max_J;
            EJ = @(J) (rot.B0*J.*(J+1))*1e2*h*c; % in the unit "J"
            Z = sum(gJ(vib_J).*(2*vib_J+1).*exp(-EJ(vib_J)/k/gas.temperature)); % partition function considering only the ground vibrational state
            rho = @(J) gJ(J).*exp(-EJ(J)/k/gas.temperature)/Z; % population
            
            % frequency shift of each rotational level
            rot.omega = (EJ(rot_J+2)-EJ(rot_J))/hbar*1e-12; % 2*pi*THz
    
            % Raman response prefactor
            rot.preR = gas.Ng(gas_i)*1/60/hbar*rot.gamma^2*(rho(rot_J)-rho(rot_J+2)).*(rot_J+2).*(rot_J+1)./(2*rot_J+3);
            
            gas.N2O = struct('R', rot);
        case 'CO2'
            % T2
            rot.T2 = 1e6/(pi*3570*eta(gas_i)); % ps; from N2, haven't investigated yet
    
            % polarizability
            % From
            % 1. Dempsey et al., "Comparative study of optical nonlinearitie of CO2 and N2O via single-shot spatio-temporally-resolved visualization," Opt. Commun. 545, 129669 (2023)
            % 2. Truong et al., "Spectral Broadening and Pulse Compression in Molecular Gas-Filled Hollow-Core Fibers," IEEE J. Sel. Topics Quantum Electron. 30(6), 1-11 (2024)
            % I feel lazy to convert values between units, so I obtain N2O's gamma by scaling it with N2's value in the second paper
            N2_paper = 9.3;
            N2O_paper = 16.3;
            rot.gamma = 4.13*au_polarizability*(N2O_paper/N2_paper); % 4.13 is from N2's
            
            % calibration factor
            % To match with the experiments of some papers
            % (I haven't tested yet, so I put one)
            rot.polarizability_calibration = 1.17;
            % Apply the calibration
            rot.gamma = rot.gamma*rot.polarizability_calibration;
            
            % nuclear spin statistical constant
            gJ = @(J) 1; % always 1
            max_J = 100; % only 33 rotational energy levels below the 1st vibrational energy level
            
            % Energy of each rotational state7
            rot.B0 = 0.3902; % cm^(-1)
            
            rot_J = 0:max_J-2;
            vib_J = 0:max_J;
            EJ = @(J) (rot.B0*J.*(J+1))*1e2*h*c; % in the unit "J"
            Z = sum(gJ(vib_J).*(2*vib_J+1).*exp(-EJ(vib_J)/k/gas.temperature)); % partition function considering only the ground vibrational state
            rho = @(J) gJ(J).*exp(-EJ(J)/k/gas.temperature)/Z; % population
            
            % frequency shift of each rotational level
            rot.omega = (EJ(rot_J+2)-EJ(rot_J))/hbar*1e-12; % 2*pi*THz
    
            % Raman response prefactor
            rot.preR = gas.Ng(gas_i)*1/60/hbar*rot.gamma^2*(rho(rot_J)-rho(rot_J+2)).*(rot_J+2).*(rot_J+1)./(2*rot_J+3);
            
            gas.CO2 = struct('R', rot);
    end
end

end