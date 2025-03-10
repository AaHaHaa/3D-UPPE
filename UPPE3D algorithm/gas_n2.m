function n2 = gas_n2(material,wavelength,eta)
%GAS_N2 Load gas's n2 value

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
num_gas = length(material);
n2 = 0; % initialization
for gas_i = 1:num_gas
    switch material{gas_i}
        case 'H2' % m^2/(W*atm)
                  % This value is taken from Wahlstrand, et al., 
                  % "Absolute measurement of the ultrafast nonlinear electronic and rovibrational response in H2 and D2" (2015)
                  % Its value, after converted into X3, is close to the paper by Belli et al.,
                  % "Vacuum-ultraviolet to infrared supercontinuum in hydrogen-filled photonic crystal fiber" Optica (2015)
                  % with X3 = 2.206e-26 m^2/V^2 at standard conditions
            n2_i = 0.65e-23;
        case 'D2' % m^2/(W*atm)
                  % This value is taken from Wahlstrand, et al., 
                  % "Absolute measurement of the ultrafast nonlinear electronic and rovibrational response in H2 and D2" (2015)
            n2_i = 0.58e-23;
        case 'N2' % m^2/(W*atm)
                  % From Jeffrey M. Brown et al.,
                  % "Ab initio calculations of the linear and nonlinear susceptibilities of N2, O2, and air in midinfrared laser pulses"
            P_n2 = 14.63e9; % W
            lambda0_n2 = 0.3334e-6; % m
            n2_i = P_n2^(-1)./(lambda0_n2^(-2) - wavelength.^(-2));
            n2_i(wavelength < 400e-9) = P_n2^(-1)./(lambda0_n2^(-2) - (400e-9).^(-2)); % this equation is valid only above 400 nm
            n2_i(isinf(n2_i)) = 0; % avoid the singularity at lambda0_n2
        case 'O2' % m^2/(W*atm)
                  % From Jeffrey M. Brown et al.,
                  % "Ab initio calculations of the linear and nonlinear susceptibilities of N2, O2, and air in midinfrared laser pulses"
            P_n2 = 14.62e9; % W
            lambda0_n2 = 0.3360e-6; % m
            n2_i = P_n2^(-1)./(lambda0_n2^(-2) - wavelength.^(-2));
            n2_i(wavelength < 400e-9) = P_n2^(-1)./(lambda0_n2^(-2) - (400e-9).^(-2)); % this equation is valid only above 400 nm
            n2_i(isinf(n2_i)) = 0; % avoid the singularity at lambda0_n2
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
            n2_i = 0.79*n2_N2 + 0.21*n2_O2; %clear n2_N2 n2_O2
            n2_i(isinf(n2_i)) = 0; % avoid the singularity
        case 'Xe' % m^2/(W*atm)
                  % From Shu-Zee Alencious Lo, et al.,
                  % "Pulse propagation in hollow-core fiber at high-pressure regime: application to compression of tens of ?J pulses and determination of nonlinear refractive index of xenon at 1.03um" Applied Optics (2018)
            n2_i = 50.1e-24;
        case 'Ar' % m^2/(W*atm)
                  % From Carsten Bree, Ayhan Demircan, and Gunter Steinmeyer,
                  % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
            n2_i = 7.96e-24;
        case 'Ne' % m^2/(W*atm)
                  % From Carsten Bree, Ayhan Demircan, and Gunter Steinmeyer,
                  % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
            n2_i = 0.85e-24;
        case 'He' % m^2/(W*atm)
                  % From Carsten Bree, Ayhan Demircan, and Gunter Steinmeyer,
                  % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
            n2_i = 0.34e-24;
        case 'Kr' % m^2/(W*atm)
                  % From Carsten Bree, Ayhan Demircan, and Gunter Steinmeyer,
                  % "Method for Computing the Nonlinear Refractive Index via Keldysh Theory" (2010)
            n2_i = 18.9e-24;
        case 'CH4'
            n2_i = 3.118e-23; % m^2/(W*atm)
        case 'N2O' % m^2/(W*atm)
                   % From Dempsey, et al.,
                   % "Comparative study of optical nonlinearities of CO2 and N2O via single-shot spatio-temporally-resolved visualization" (2023)
            n2_i = 21e-24;
        case 'CO2' % m^2/(W*atm)
                   % From Dempsey, et al.,
                   % "Comparative study of optical nonlinearities of CO2 and N2O via single-shot spatio-temporally-resolved visualization" (2023)
            n2_i = 13e-24;
    end
    n2 = n2 + n2_i*eta(gas_i);
end

end

