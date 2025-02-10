% This code reads the pressure-induced absorption data of O2 from the
% following paper:
% Crawford et al., "Infra-Red Absorption of Oxygen and Nitrogen Induced by Intermolecular Forces," Phys. Rev. 75, 1607 (1949)
%
% Although Raman should be IR inactive, molecular collision can distort
% electron distribution resulting in IR-active Raman absorption.

clearvars; close all;

data = readmatrix('O2 vib Raman absorption from collision.csv');

[wavenumber,idx] = sort(data(:,1)); % cm^(-1)
absorption = smooth(data(idx,2),3); % cm^(-1)*amagat^(-2)

%% Extend the wavenumber window
% Left
extended_wavenumber = (1:1352)'; % cm^(-1)
% B.C.: Slope at the edge
absorption_slope_left = (absorption(4)-absorption(1))/(wavenumber(4)-wavenumber(1));
% Compute extrapolations
extended_absorption = exponential_decay(wavenumber(1),absorption(1),absorption_slope_left,extended_wavenumber);
% extended_part
wavenumber = [extended_wavenumber;wavenumber];
absorption = [extended_absorption;absorption];

% Right
extended_wavenumber = (1757:1e4)'; % cm^(-1)
% B.C.: Slope at the edge
absorption_slope_right = (absorption(end)-absorption(end-3))/(wavenumber(end)-wavenumber(end-3));
% Compute extrapolations
extended_absorption = exponential_decay(wavenumber(end),absorption(end),absorption_slope_right,extended_wavenumber);
% extended_part
wavenumber = [wavenumber;extended_wavenumber];
absorption = [absorption;extended_absorption];

%%

save('O2_vib_Raman_absorption_from_collision.mat','wavenumber','absorption');

figure;
plot(wavenumber,absorption,'linewidth',2);
xlabel('Wavenumber (cm^{-1})');
ylabel('\alpha (cm^{-1}amagat^{-2})');
set(gca,'fontsize',20);