function absorption = read_absorption(material,lambda,gas_density)
%READ_ABSORPTION It reads the pressure-induced Raman-IR absorption.
% Although Raman should be IR inactive, molecular collision can distort
% electron distribution resulting in IR-active Raman absorption.
%
% Input:
%   material: 'H2','D2','N2','O2','air'
%   lambda: wavelength (m)
%   gas_density: gas density in amagat

% "lambda" needs to be a column vector.
if size(lambda,1)==1
    lambda = lambda.';
end
% "lambda" needs to be in the order of small to large because of extrapolation below
% If not, sort it here and the result will be reverse back to the unsorted original one after computations.
if length(lambda)>1
    lambda_sort = true;
    [lambda,sort_idx] = sort(lambda);
    [~,reverse_sort_idx] = sort(sort_idx);
else
    lambda_sort = false;
end

%% Reading data from the specified file
switch material
    case 'H2'
        dataArray_rot = load('H2_rot_Raman_absorption_from_collision.mat');
        dataArray_vib = load('H2_vib_Raman_absorption_from_collision.mat');
        
        wavenumber_raw = [dataArray_rot.wavenumber;dataArray_vib.wavenumber]; % cm^(-1)
        absorption_raw = [dataArray_rot.absorption;dataArray_vib.absorption]; % cm^(-1)*amagat^(-2)
        
        absorption_raw = absorption_raw*gas_density^2;
    case 'D2'
        dataArray = load('D2_Raman_absorption_from_collision.mat');
        
        wavenumber_raw = [dataArray.wavenumber]; % cm^(-1)
        absorption_raw = [dataArray.absorption]; % cm^(-1)*amagat^(-2)
        
        absorption_raw = absorption_raw*gas_density^2;
    case 'N2'
        dataArray_vib = load('N2_vib_Raman_absorption_from_collision.mat');
        
        wavenumber_raw = dataArray_vib.wavenumber; % cm^(-1)
        absorption_raw = dataArray_vib.absorption; % norm.
        
        a1 = 3.83e-4; % cm^(-2)*amagaat^(-2)
        a2 = -2.9e-7; % cm^(-2)*amagaat^(-3)
        integrated_absorption = a1*gas_density^2 + a2*gas_density^3;
        
        absorption_raw = absorption_raw*(integrated_absorption/trapz(wavenumber_raw,absorption_raw)); % cm^(-1)*amagat^(-2)
    case 'O2'
        dataArray_vib = load('O2_vib_Raman_absorption_from_collision.mat');
        
        wavenumber_raw = dataArray_vib.wavenumber; % cm^(-1)
        absorption_raw = dataArray_vib.absorption; % cm^(-1)*amagat^(-2)
        
        absorption_raw = absorption_raw*gas_density^2;
    case 'air'
        % N2
        dataArray_vib = load('N2_vib_Raman_absorption_from_collision.mat');
        
        wavenumber_raw_N2 = dataArray_vib.wavenumber; % cm^(-1)
        absorption_raw_N2 = dataArray_vib.absorption; % norm.
        
        a1 = 3.83e-4; % cm^(-2)*amagaat^(-2)
        a2 = -2.9e-7; % cm^(-2)*amagaat^(-3)
        integrated_absorption = a1*gas_density^2 + a2*gas_density^3;
        
        absorption_raw_N2 = absorption_raw_N2*(integrated_absorption/trapz(wavenumber_raw_N2,absorption_raw_N2)); % cm^(-1)*amagat^(-2)
        
        % O2
        dataArray_vib = load('O2_vib_Raman_absorption_from_collision.mat');
        
        wavenumber_raw_O2 = dataArray_vib.wavenumber; % cm^(-1)
        absorption_raw_O2 = dataArray_vib.absorption; % cm^(-1)*amagat^(-2)
        
        absorption_raw_O2 = absorption_raw_O2*gas_density^2;
        
        wavenumber_raw = sort(unique([wavenumber_raw_O2;wavenumber_raw_N2]));
        absorption_raw = interp1(wavenumber_raw_N2,absorption_raw_N2,wavenumber_raw,'pchip')*0.8 + ...
                         interp1(wavenumber_raw_O2,absorption_raw_O2,wavenumber_raw,'pchip')*0.2;
end
if wavenumber_raw(1) == 0
    wavenumber_raw(1) = 1; % cm^(-1)
end

%% Defining the raw data (from the file)
lambda_raw = 1e-2./wavenumber_raw; % m
absorption_raw = absorption_raw*1e2; % m^(-1)*amagat^(-2)

lambda_raw = flipud(lambda_raw);
absorption_raw = flipud(absorption_raw);

%% Interpolating data for our wavelength grid
lambda_raw_min = lambda_raw(1);
lambda_raw_max = lambda_raw(end);
lambda_inside_lambda_raw = lambda(lambda<=lambda_raw_max & lambda>=lambda_raw_min); % interpolation
lambda_outside_lambda_raw_left  = lambda(lambda<lambda_raw_min); % extrapolation
lambda_outside_lambda_raw_right = lambda(lambda>lambda_raw_max); % extrapolation

% interpolation
absorption = [zeros(size(lambda_outside_lambda_raw_left)); ...
                interp1(lambda_raw,absorption_raw, lambda_inside_lambda_raw, 'pchip'); ...
              zeros(size(lambda_outside_lambda_raw_right))];

% extrapolation with an exponential decay
% The left part
if ~isempty(lambda_outside_lambda_raw_left)
    % B.C.: Slope at the edge
    absorption_raw_slope_left = (absorption_raw(2)-absorption_raw(1))/(lambda_raw(2)-lambda_raw(1));
    % Compute extrapolations
    absorption(1:length(lambda_outside_lambda_raw_left)) = exponential_decay(lambda_raw(1),absorption_raw(1),max(0,absorption_raw_slope_left),lambda_outside_lambda_raw_left);
end
% The right part
if ~isempty(lambda_outside_lambda_raw_right)
    % B.C.: Slope at the edge
    absorption_raw_slope_right = (absorption_raw(end)-absorption_raw(end-1))/(lambda_raw(end)-lambda_raw(end-1));
    % Compute extrapolations
    absorption(end-length(lambda_outside_lambda_raw_right)+1:end) = exponential_decay(lambda_raw(end),absorption_raw(end),min(0,absorption_raw_slope_right),lambda_outside_lambda_raw_right);
end

% Reverse back the original order of "lambda"
if lambda_sort
    absorption = absorption(reverse_sort_idx);
end

end