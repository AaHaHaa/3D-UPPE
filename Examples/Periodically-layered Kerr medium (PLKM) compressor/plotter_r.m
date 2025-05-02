function fig = plotter_r(fig,...
                         A,...
                         z,MFD,...
                         r,lambda)

% Make them column vectors
z = z(:);
MFD = MFD(:);

plot_wavelength_lim = [950,1100];
plot_r_lim = [0,220];

spectrum = abs(fftshift(ifft(A),1)).^2./lambda.^2; % "./lambda" is to make it into the wavelength domain (not with the correct unit but the correct relative strength; we'll normalize it later)

if isempty(fig)
    fig = figure;
else
    figure(fig); % reuse the previous figure handle; plot in the same figure
end
subplot(2,2,1)
pcolor(lambda,r*1e6,spectrum(:,:,1)');
shading interp; colormap(jet);
xlabel('Wavelength (nm)');
ylabel('r (\mum)');
xlim(plot_wavelength_lim);
ylim(plot_r_lim);
title('Beam input');

subplot(2,2,2)
pcolor(lambda,r*1e6,spectrum(:,:,end)');
shading interp; colormap(jet);
xlabel('Wavelength (nm)');
ylabel('r (\mum)');
xlim(plot_wavelength_lim);
ylim(plot_r_lim);
title('Beam output');

subplot(2,2,3)
plot(z*1e2,MFD*1e3,'Color','k','linewidth',2);
xlabel('z (cm)');
title('Beam size evolution')

subplot(2,2,4)
r_idx = r*1e3 < MFD(end)/2;
center_spectrum = spectrum(:,1,end);
avg_spectrum = trapz(r(r_idx),spectrum(:,r_idx,end).*r(r_idx),2);
avg_spectrum = avg_spectrum/max(avg_spectrum); % normalized
center_spectrum = center_spectrum/max(center_spectrum); % normalized
plot(lambda,avg_spectrum,'Color','b','linewidth',2);
hold on;
plot(lambda,center_spectrum,'Color','r','linewidth',2);
hold off;
xlabel('Wavelength (nm)');
xlim(plot_wavelength_lim);
legend('Avg spectrum','Center spectrum');
title('Spectral domain')

drawnow; % plot it during running the code

end