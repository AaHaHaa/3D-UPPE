clearvars; close all;

f = linspace(0.01,200,10000)'; % THz
lambda = 299792458./f*1e-12; % m
gas_density = 1; % amagat

absorption_H2 = read_absorption('H2',lambda,gas_density);
absorption_N2 = read_absorption('N2',lambda,gas_density);
absorption_O2 = read_absorption('O2',lambda,gas_density);

figure;
tl = tiledlayout(1,1);
ax1 = axes(tl);
h = plot(ax1,1./lambda/100,[absorption_H2,absorption_N2,absorption_O2]*1e4,'linewidth',2);
l = legend('H_2','N_2','O_2'); legend('boxoff');
set(h(1),'Color','b'); set(h(2),'Color','r'); set(h(3),'Color','k');
xlim([0,6500]); set(gca,'fontsize',20);
xlabel('Wavenumber (_{}cm^{-1})'); box off;
xtick = get(ax1,'XTick'); xlimits = get(ax1,'XLim'); ylimits = get(ax1,'YLim');
xtick_ratio = xtick/xlimits(2);
lambdatick = [50,10,5,3,2,1.6]; lambdatick_ratio = 1e4./lambdatick/xlimits(2);
ylabel('Absorption (10^{-4}m^{-1}amagat^{-2})');
ax2 = axes(tl,'XAxisLocation','top','Box','off','color','none','YAxisLocation','right');
set(ax2,'YLim',ylimits,'XTickLabel',arrayfun(@(x)num2str(x),lambdatick,'UniformOutput',false),'XTick',lambdatick_ratio,'YTick',[]);
xlabel('Wavelength (µm)');
set(gca,'fontsize',20);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),560,470]);
print('pressure induced absorption (H2 N2 O2).eps','-depsc');