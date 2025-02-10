function [D_op_upsampling,W_op_upsampling,loss_op_upsampling,...
          prefactor,...
          upsampling_zeros] = upsampling(sim,...
                                         Nt,gas_Nt,...
                                         Nr,Np,...
                                         D_op,W_op,loss_op,...
                                         prefactor)
%UPSAMPLING Upsampling to avoid frequency aliasing
%
% Aliasing mostly comes from Raman shift. However, when the spectrum
% becomes broad, aliasing can come from Kerr effect as well due to
% four-wave mixing.
% Instead of applying upsampling to the Raman computation only, it's
% important to apply it to Kerr as well, especially when running
% supercontinuum generation or when your frequency window isn't large
% enough.
if sim.gpu_yes
    upsampling_zeros = complex(zeros(gas_Nt-Nt, Nr, 1, Np, 'gpuArray'));
else
    upsampling_zeros = complex(zeros(gas_Nt-Nt, Nr, 1, Np));
end

n = ceil(Nt/2);
% Zero-paddding in frequency domain for upsampling temporally
% Zeros are added at the low- and high-frequency side which is the center
% of the array after discrete Fourier transform.
D_op_upsampling    = cat(1,D_op(1:n,:,:,:),   zeros([gas_Nt-Nt, size(D_op,2:4)]),   D_op(n+1:end,:,:,:));
W_op_upsampling    = cat(1,W_op(1:n,:,:,:),   zeros([gas_Nt-Nt, size(W_op,2:4)]),   W_op(n+1:end,:,:,:));
loss_op_upsampling = cat(1,loss_op(1:n,:,:,:),zeros([gas_Nt-Nt, size(loss_op,2:4)]),loss_op(n+1:end,:,:,:));

prefactor{1} = cat(1,prefactor{1}(1:n,:,:,:),zeros([gas_Nt-Nt, size(prefactor{1},2:4)]),prefactor{1}(n+1:end,:,:,:));
prefactor{2} = cat(1,prefactor{2}(1:n,:,:,:),zeros([gas_Nt-Nt, size(prefactor{2},2:4)]),prefactor{2}(n+1:end,:,:,:));
prefactor{3} = cat(1,prefactor{3}(1:n,:,:,:),zeros([gas_Nt-Nt, size(prefactor{3},2:4)]),prefactor{3}(n+1:end,:,:,:));

if sim.photoionization_model
    prefactor{4} = cat(1,prefactor{4}(1:n,:,:,:),zeros([gas_Nt-Nt, size(prefactor{4},2:4)]),prefactor{4}(n+1:end,:,:,:));
    prefactor{5} = cat(1,prefactor{5}(1:n,:,:,:),zeros([gas_Nt-Nt, size(prefactor{5},2:4)]),prefactor{5}(n+1:end,:,:,:));
end

end

