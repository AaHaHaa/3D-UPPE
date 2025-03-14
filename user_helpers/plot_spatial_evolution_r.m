function varargout = plot_spatial_evolution_r(r,z,A2,varargin)
%PLOT_SPATIAL_EVOLUTION_R Summary of this function goes here
%   Detailed explanation goes here

r_max = max(r); % m

Nr = length(r);
Nx = Nr*2 - 1;
x = linspace(-r_max,r_max,Nx)'; % m

A2 = interp1(r.',A2,abs(x),'linear',0);

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

