function Frame = animator_r(Frame,...
                            A,...
                            z,MFD,start_idx,...
                            Nt,dt,r,lambda,...
                            plate_z)

% Make them column vectors
z = z(:);
MFD = MFD(:);

c = 299792.458; % nm/ps
factor_correct_unit = (Nt*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                     % "/1e3" is to make pJ into nJ
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

plot_wavelength_lim = [950,1100];
plot_MFD_lim = [400,620];
for j = 1:size(A,3)-1
    if exist('fig','var')
        figure(fig);
    else
        fig = figure;
    end

    spectrum = abs(fftshift(ifft(A(:,:,j+1),[],1),1)).^2*factor_correct_unit.*factor; % nJ/nm/m^2

    subplot(2,2,1);
    pcolor(lambda,r*1e6,spectrum');
    shading interp; colormap(jet);
    xlabel('Wavelength (nm)');
    ylabel('r (\mum)');
    xlim(plot_wavelength_lim);
    ylim([0,plot_MFD_lim(2)/2]);

    subplot(2,2,2);
    plot(plate_z(1)*1e2*[1;1],plot_MFD_lim,'Color','b','linewidth',1.5);
    hold on;
    for i = 2:length(plate_z)
        switch mod(i-1,2)
            case 1
                plot(plate_z(i)*1e2*[1;1],plot_MFD_lim,'Color','m','linewidth',0.5);
            case 0
                plot(plate_z(i)*1e2*[1;1],plot_MFD_lim,'Color','b','linewidth',1.5);
        end
        hold on;
    end
    plot(z(1:start_idx+j)*1e2,MFD(1:start_idx+j)*1e3,'Color','k','linewidth',2);
    hold off;
    xlabel('z (cm)');
    ylabel('MFD (\mum)');
    xlim([0,plate_z(end)*1e2]); % cm
    ylim(plot_MFD_lim);

    subplot(2,2,[3,4]);
    r_idx = r*1e3 < MFD(end)/2;
    center_spectrum = spectrum(:,1);
    avg_spectrum = trapz(r(r_idx),spectrum(:,r_idx).*r(r_idx),2);
    avg_spectrum = avg_spectrum/max(avg_spectrum); % normalized
    center_spectrum = center_spectrum/max(center_spectrum); % normalized
    plot(lambda,avg_spectrum,'Color','b','linewidth',2);
    hold on;
    plot(lambda,center_spectrum,'Color','r','linewidth',2);
    hold off;
    xlabel('Wavelength (nm)');
    ylabel('PSD (norm.)');
    xlim(plot_wavelength_lim);
    ylim([0,1]);
    legend('Avg spectrum','Center spectrum (norm.)');

    set(fig,'Color',[1,1,1]);

    Frame(j) = getframe(fig);
end
close(fig);

end