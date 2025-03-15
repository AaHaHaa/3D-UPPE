close all;

x = linspace(0,1,10000)';

[K,E] = ellipke(x);

figure;
h = plot(x,[K,E],'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r');
set(gca,'fontsize',20);
xlabel('x');
l = legend('K(x)','E(x)'); set(l,'fontsize',20,'location','northwest');
print(gcf,'KE.pdf','-dpdf');