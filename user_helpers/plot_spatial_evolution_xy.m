function varargout = plot_spatial_evolution_xy(x,z,A2,varargin)
%PLOT_SPATIAL_EVOLUTION_XY Summary of this function goes here
%   Detailed explanation goes here

if isempty(varargin)
    fig = figure;
else
    fig = varargin{1};
    figure(varargin{1});
end
pcolor(z,x,A2);
shading interp; colormap(jet);
xlabel('Propagation (m)');
ylabel('X (m)');
set(gca,'fontsize',20);
ax = gca;

varargout = {fig,ax};

end

