% This code reads the pressure-induced absorption data of N2 from the
% following paper:
% S. Paddi Reddy and C. W. Cho, "INDUCED INFRARED ABSORPTION OF NITROGEN AND NITROGEN – FOREIGN GAS MIXTURES," Canadian Journal of Physics 43, 2331-2343 (1965)
%
% Although Raman should be IR inactive, molecular collision can distort
% electron distribution resulting in IR-active Raman absorption.
%
% The data is normalized, so we need re-normalize it to the actual values
% a1 = 3.83e-4; % cm^(-2)*amagaat^(-2)
% a2 = -2.9e-7; % cm^(-2)*amagaat^(-3)
% integrated_absorption = a1*eta^2 + a2*eta^3, eta is the gas density in amagat
% integrated absorption is integral(wavenumber,absorption)

clearvars; close all;

data = readmatrix('N2 vib Raman absorption from collision.csv');

[wavenumber,idx] = sort(data(:,1)); % cm^(-1)
absorption = smooth(data(idx,2),3); % cm^(-1)*amagat^(-2)

%% Extend the wavenumber window
% Left
extended_wavenumber = (1:2060)'; % cm^(-1)
% B.C.: Slope at the edge
absorption_slope_left = (absorption(4)-absorption(1))/(wavenumber(4)-wavenumber(1));
% Compute extrapolations
extended_absorption = exponential_decay(wavenumber(1),absorption(1),absorption_slope_left,extended_wavenumber);
% extended_part
wavenumber = [extended_wavenumber;wavenumber];
absorption = [extended_absorption;absorption];

% Right
extended_wavenumber = (2705:1e4)'; % cm^(-1)
% B.C.: Slope at the edge
absorption_slope_right = (absorption(end)-absorption(end-3))/(wavenumber(end)-wavenumber(end-3));
% Compute extrapolations
extended_absorption = exponential_decay(wavenumber(end),absorption(end),absorption_slope_right,extended_wavenumber);
% extended_part
wavenumber = [wavenumber;extended_wavenumber];
absorption = [absorption;extended_absorption];

%%

save('N2_vib_Raman_absorption_from_collision.mat','wavenumber','absorption');

figure;
plot(wavenumber,absorption,'linewidth',2);
xlabel('Wavenumber (cm^{-1})');
ylabel('\alpha (norm.)');
set(gca,'fontsize',20);

gas_density = 1; % plot with 1 amagat
absorption = read_absorption('N2',1e-2./wavenumber,gas_density)*1e-2; % cm^(-1)
figure;
plot(wavenumber,absorption,'linewidth',2);
xlabel('Wavenumber (cm^{-1})');
ylabel('\alpha (cm^{-1}amagat^{-2})');
set(gca,'fontsize',20);

%% Verification
fprintf('integrated absorption (this data):    %3.2e\n',trapz(wavenumber,absorption));

a1 = 3.83e-4; % cm^(-2)*amagaat^(-2)
a2 = -2.9e-7; % cm^(-2)*amagaat^(-3)
integrated_absorption = a1*gas_density^2 + a2*gas_density^3;
fprintf('integrated absorption (in the paper): %3.2e\n',integrated_absorption);