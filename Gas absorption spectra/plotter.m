clearvars; close all;

f = linspace(0.01,200,10000)'; % THz
lambda = 299792458./f*1e-12; % m
gas_density = 1; % amagat
absorption = read_absorption('O2',lambda,gas_density);

figure;
tl = tiledlayout(1,1);
ax1 = axes(tl);
plot(ax1,1./lambda/100,absorption*1e4,'linewidth',2,'Color','b');
xlim([0,6500]); set(gca,'fontsize',20);
xlabel('Wavenumber (_{}cm^{-1})'); box off;
xtick = get(ax1,'XTick'); xlimits = get(ax1,'XLim'); ylimits = get(ax1,'YLim');
xtick_ratio = xtick/xlimits(2);
lambdatick = [50,10,5,3,2,1.6]; lambdatick_ratio = 1e4./lambdatick/xlimits(2);
ylabel('Absorption (10^{-4}m^{-1}amagat^{-2})');
ax2 = axes(tl,'Position',get(ax1,'Position'),'XAxisLocation','top','Box','off','color','none','YAxisLocation','right');
set(ax2,'YLim',ylimits,'XTickLabel',arrayfun(@(x)num2str(x),lambdatick,'UniformOutput',false),'XTick',lambdatick_ratio,'YTick',[]);
xlabel('Wavelength (µm)');
set(gca,'fontsize',20);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),560,470]);
%print('pressure induced absorption.eps','-depsc');