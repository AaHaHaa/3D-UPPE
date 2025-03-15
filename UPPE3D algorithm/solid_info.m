function [fiber,sim,solid] = solid_info(fiber,sim,solid)
%SOLID_INFO It loads parameters related to solids

%% Error check
if sim.gpu_yes
    try
        gpuDevice;
    catch
        error('solid_info_constant_pressure:GPUError',...
              'No GPU is detected. Please set "sim.gpu_yes=false".');
    end
end

%% Nonlinear coefficient
fiber.n2 = solid_n2(solid.material{1});

%% Ionization potential
if sim.photoionization_model
    solid = photoionization_parameters(solid);
end

%%
% There will be some situations where we want to modify, for example, the 
% Raman response functions, so solid_info() will be called before pulse
% propagation. In this case, we don't want the propagation function to call
% solid_info() again and rewrite the desired parameters set beforehand.
solid.info_called = true; % info() has been called

end