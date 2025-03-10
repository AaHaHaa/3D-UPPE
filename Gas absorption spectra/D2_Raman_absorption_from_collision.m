% This code reads the pressure-induced absorption data of D2 from the
% following paper:
% Abel and Frommhold, "Note: Collision-induced infrared absorption by gaseous deuterium," J. Chem. Phys. 133, 146101(2010)
%
% I extracted the one at 200 K.
%
% Although Raman should be IR inactive, molecular collision can distort
% electron distribution resulting in IR-active Raman absorption.

clearvars; close all;

data = readmatrix('D2 Raman absorption from collision.csv');

[wavenumber,idx] = sort(data(:,1)); % cm^(-1)
absorption = data(idx,2); % cm^(-1)*amagat^(-2)

%%

save('D2_Raman_absorption_from_collision.mat','wavenumber','absorption');

figure;
plot(wavenumber,absorption,'linewidth',2);
xlabel('Wavenumber (cm^{-1})');
ylabel('\alpha (cm^{-1}amagat^{-2})');
set(gca,'fontsize',20);