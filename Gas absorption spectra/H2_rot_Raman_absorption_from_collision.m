% This code reads the pressure-induced absorption data of H2 from the
% following paper:
% Kiss et al, "THE PRESSURE-INDUCED ROTATIONAL ABSORPTION SPECTRUM OF HYDROGEN: I. A STUDY OF THE ABSORPTION INTENSITIES," Canadian Journal of Physics 37, 362-376 (1959) 
%
% Although Raman should be IR inactive, molecular collision can distort
% electron distribution resulting in IR-active Raman absorption.

clearvars; close all;

data = readmatrix('H2 rot Raman absorption from collision.csv');

data(:,2) = data(:,2)/max(data(:,2))*4.19e-6; % normalized to the max of 4.19e-6 (actual number)

[wavenumber,idx] = sort(data(:,1)); % cm^(-1)
absorption = data(idx,2); % cm^(-1)*amagat^(-2)

%% Extend the wavenumber window
% Left
extended_wavenumber = (1:290)'; % cm^(-1)
% B.C.: Slope at the edge
absorption_slope_left = (absorption(4)-absorption(1))/(wavenumber(4)-wavenumber(1));
% Compute extrapolations
extended_absorption = exponential_decay(wavenumber(1),absorption(1),absorption_slope_left,extended_wavenumber);
% extended_part
wavenumber = [extended_wavenumber;wavenumber];
absorption = [extended_absorption;absorption];

% Right
extended_wavenumber = (1400:2000)'; % cm^(-1)
% B.C.: Slope at the edge
absorption_slope_right = (absorption(end)-absorption(end-3))/(wavenumber(end)-wavenumber(end-3));
% Compute extrapolations
extended_absorption = exponential_decay(wavenumber(end),absorption(end),absorption_slope_right,extended_wavenumber);
% extended_part
wavenumber = [wavenumber;extended_wavenumber];
absorption = [absorption;extended_absorption];

%%

save('H2_rot_Raman_absorption_from_collision.mat','wavenumber','absorption');

figure;
plot(wavenumber,absorption,'linewidth',2);
xlabel('Wavenumber (cm^{-1})');
ylabel('\alpha (cm^{-1}amagat^{-2})');
set(gca,'fontsize',20);

%% Verification
fprintf('integrated absorption @1 amagat (this data):    %3.2e\n',trapz(wavenumber,absorption));
fprintf('integrated absorption @1 amagat (in the paper): %3.2e\n',2.06e-3);