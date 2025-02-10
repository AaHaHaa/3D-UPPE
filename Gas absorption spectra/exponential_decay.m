function y = exponential_decay(x_edge,A,B,x)
%EXPONENTIAL_DECAY
%   It extrapolates with the function exp((+-)a*(x-b)).
%
%   B.C.:
%       continuity of points: exp((+-)a*(x_edge-b))       = A
%       continuity of slope:  (+-)a*exp((+-)a*(x_edge-b)) = B
%   
%   solution:
%       (+-)a = B/A
%       b = x_edge - ((+-)log(A))/a = x_edge - A*log(A)/B

if A == 0
    y = zeros(size(x));
elseif B == 0
    % choose somewhere near the center of spectrum as "b" and find "a"
    b = 1e-6; % wavelength: 1um
    pma = log(A)/(x_edge-b);
    
    % extrapolation
    y = exp(pma*(x-b));
else
    % coefficients of exponential decay function
    pma = B/A;
    b = x_edge - A*log(A)/B;

    % extrapolation
    y = exp(pma*(x-b));
end

y = y.*exp(-abs(x-x_edge)/300); % make it decay within 300 cm^(-1)

end