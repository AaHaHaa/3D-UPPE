function foutput = UPPE3D_propagate(fiber, initial_condition, sim, gas)
%UPPE3D_PROPAGATE Propagate an initial full-field pulse through an arbitrary 
% distance of a nonlinear medium, such as an optical fiber
%   This is a caller function, calling
%   UPPE_propagate_xy(),
%   UPPE_propagate_r(), or
%   UPPE_propagate_r_gas()
%   based on whether the field is radially symmetric and whether the medium is gas or not.

%%
if ispc
    sep_char = '\';
else % unix
    sep_char = '/';
end

% Load the folder
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep_char);
upper_folder = current_path(1:sep_pos(end-1));

sim.cuda_dir_path = [upper_folder 'cuda'];

%% Determine whether to use the radially-symmetric scheme
rxy_str = 'xy';
if isfield(initial_condition,'r')
    rxy_str = 'r';
end

gas_str = '';
if exist('gas','var')
     gas_str = '_gas';

     if isequal(rxy_str,'xy')
         error('UPPE3D_propagate:simulationSchemeError',...
               'Gas-filled simulation supports only the radially-symmetric scheme.');
     end
else
    gas = [];
end

UPPE3D_propgation_func = str2func(['UPPE3D_propagate_', rxy_str,gas_str]);

%% Run the pulse propagation
foutput = UPPE3D_propgation_func(fiber, initial_condition, sim, gas);

end