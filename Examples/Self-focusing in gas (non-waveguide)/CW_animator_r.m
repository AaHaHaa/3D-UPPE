function Frame = CW_animator_r(A,...
                               z,MFD,...
                               Nt,dt,r,lambda,...
                               L0)

% Make them column vectors
z = z(:);
MFD = MFD(:);

c = 299792.458; % nm/ps
factor_correct_unit = (Nt*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                     % "/1e3" is to make pJ into nJ
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

Frame(size(A,3)) = struct('cdata',[],'colormap',[]);
for j = 1:size(A,3)
    if exist('fig','var')
        figure(fig);
    else
        fig = figure;
    end

    subplot(1,2,1);
    A_spatial = A(floor(Nt/2)+1,:,j);
    plot(r*1e6,abs(A_spatial.').^2/max(abs(A_spatial).^2),'linewidth',2,'Color','b');
    xlabel('r (\mum)');
    xlim([0,1]*1e2);

    subplot(1,2,2);
    plot(z(1:j)*1e2,MFD(1:j),'Color','k','linewidth',2);
    xlabel('z (cm)');
    ylabel('MFD (\mum)');
    xlim([0,L0*1e2]); % cm
    ylim([0,100]); % um

    set(fig,'Color',[1,1,1]);

    Frame(j) = getframe(fig);
end
close(fig);

end