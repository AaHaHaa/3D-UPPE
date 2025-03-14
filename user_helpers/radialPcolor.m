function varargout = radialPcolor(r,A2,varargin)
%RADIALPCOLOR Plot the radially-symmetric field
%
% r: radial sampled positions (m)
% A2: the field intensity profile with size (1,Nr)

%%
r_max = max(r); % m

Nr = length(r);
Nx = Nr*2 - 1;
x = linspace(-r_max,r_max,Nx); % m
[xx,yy] = meshgrid(x,x);

rr = sqrt(xx.^2+yy.^2);

A2 = reshape(interp1(r.',A2.',rr(:),'linear',0),Nx,Nx);

if isempty(varargin)
    fig = figure;
else
    fig = varargin{1};
    figure(varargin{1});
end
pcolor(x,x,A2);
shading interp; colormap(jet);
ax = gca;

varargout = {fig,ax};

end