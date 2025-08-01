# 3D-UPPE
This is the shared package to simulate, with MATLAB, pulse propagation in bulk crystal/gas/free space with 3D-UPPE.

It is useful for simulating solid-based or gas-filled multipass cell or multiplate compressor, etc.

## Capabilities:
1. It solves the pulse propagation with a nested [RK4IP](http://www.sciencedirect.com/science/article/pii/S0010465512004262) (Runge-Kutta under the interaction picture). Please find details in the 3D-UPPE's readme.
2. Adaptive step-size control is implemented (for the nested RK4IP).
3. Support broadband scenarios by having $\beta(\omega)$.
4. Support both scalar and polarized scenarios, controlled with `sim.scalar=true/false`.
7. Support noise-seeded processes, such as spontaneous Raman scattering, with [the newly-developed noise model](https://doi.org/10.48550/arXiv.2410.20567).
8. Efficient GPU computation (with Nvidia CUDA) is implemented. It is controlled by `sim.gpu_yes=true/false`.
9. Support radially-symmetric scheme with the Hankel transform for efficient modeling.
10. Support both solid and gas environments.
11. For gases, it supports both noble and Raman-active gases with the newly-developed vector Raman model [[1]](#references-our-papers).
12. Support photoionization in both solids and gases with the Perelomov-Popov-Terent'ev (PPT) model.
13. Support both pulsed and CW ($N_t=1$) cases. Full nonlinearity is supported in CW to fast simulate some phenomena, such as high-average-power self-focusing.

## Fourier and Hankel transforms
Since I've seen many misuse of Fourier Transform, I wrote [this tutorial](https://doi.org/10.48550/arXiv.2412.20698). Please take a look. Briefly speaking for one misuse, it's necessary to use MATLAB's `ifft` for Fourier Transform into the spectral domain.  
In addition, I have improved and implemented a new numerical Hankel transform scheme based on FHATHA, which might be publishable to a small journal (but I'm lazy). I put it in this arXiv tutorial as well. You can take a look if interested.

## How to activate CUDA for GPU computing in MATLAB:
Typically MATLAB deals with this, but there are still come steps to follow before CUDA can really be used, especially when compiling .cu files to generate .ptx files. Below I show only steps for Windows. For linux, please search for their specific steps. I've never used Mac, so I cannot comment anything on this; some functions need to be revised for extended capabilities for Mac as far as I know.<br>
1. Install [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit)
2. Install [Visual Studio Community](https://visualstudio.microsoft.com/vs/community/). Only **Desktop development with C++** is required. If it later says that it needs to install some other components due to the dependency issues, also install them.
![VS installation screenshot](Readme_images/VS_install.png)
3. Add required path of Visual Studio to computer's environmental PATH; otherwise, MATLAB, during compiling .cu files, will say "cl.exe" cannot be found.
![add PATH screenshot](Readme_images/add_PATH.png)
4. Restart the computer if something is wrong. Connections between MATLAB and CUDA or Visual Studio requires restarting to be effective.
> [!WARNING]
> MATLAB supports only a certain version of CUDA and GPUs ([support list](https://www.mathworks.com/help/releases/R2021b/parallel-computing/gpu-support-by-release.html)). CUDA or GPU that is too old just isn't supported.

## References (our papers):
1. [Raman scattering](https://doi.org/10.1063/5.0189749)
2. [Noise modeling](https://doi.org/10.48550/arXiv.2410.20567)

## Demonstrations:
- **Gas-filled Multipass cell**  
A multipass cell acts as a discrete waveguide that is commonly used for pulse compression.  
Below is an Ar-filled multipass cell that compresses a 210-fs pulse to 35 fs (dechirped pulse shown on the right).
Source: "3D-UPPE/Examples/Multipass cell/Gas-filled MPC/Ar"  
<img src="Readme_images/MPC_r.gif" width=45%><img src="Readme_images/MPC_dechirped.jpg" width=45%>

- **Periodically-layered Kerr medium**  
[Periodically-layered medium](https://doi.org/10.1364/OL.539381) in air can be a waveguide in nonlinear conditions.  
It acts as a discrete "nonlinear" waveguide with interleaving media of low (e.g., air) and high (e.g., thin glass) nonlinear refractive indices. Thin glass induces nonlinear self-focusing and air introduces diffraction. This artificially-contructed waveguide nonlinearly broadens the pulse, introducing self-phase modulation that can be compensated with a dechirper. This results in a temporally-compressed pulse. Typical compression factor is around 5.  
Source: "3D-UPPE/Examples/Periodically-layered Kerr medium (PLKM) compressor"  
<img src="Readme_images/PLKM_r.gif" width=45%><img src="Readme_images/PLKM_dechirped.jpg" width=45%>

- **Self focusing (in silica)**  
Pulse with high peak power experiences self-focusing in a Kerr medium with a positive nonlinear refractive index. The medium effectively acts as a lens, reducing the propagating beam size.  
Source: "3D-UPPE/Examples/Tutorial/2. pulsed/Self-focusing (non-waveguide)"  
<img src="Readme_images/self-focusing_r.gif" width=45%>

- **Self focusing with/without photoionization (in N2)**  
Pulse with high peak power experiences self-focusing in a gas medium as in solids. Here in this case, the peak power is made so high that it ionizes the gas, which defocuses the beam.  
Below I show examples without (left column) and with (right column) the photoionization contribution. Without it, the beam only self-focuses due to electronic nonlinearity. In addition, photoionization-induced blueshift can be clearly seen.  
Source: "3D-UPPE/Examples/Self-focusing in gas (non-waveguide)"  
<img src="Readme_images/noPhotoionization_self-focusing.jpg" width=45%><img src="Readme_images/withPhotoionization_self-focusing.jpg" width=45%>  
<img src="Readme_images/pulsed_self-focusing_noPhotoionization.gif" width=45%><img src="Readme_images/pulsed_self-focusing.gif" width=45%>

## Notes:
There is a `readme.pdf` in the `Documentations/` folder of **3D-UPPE**. Please find details of how to use this package in it. However, the fastest way to learn how to use this package is to learn from the examples in the `Examples/` folder.

I'm Yi-Hao Chen, the author of the code and from Frank Wise's group at Cornell Applied Physics.

## History:
* 2/20/2025:<br>
Ar refractive index was wrong! I fixed it.
* 3/10/2025:<br>
Added full nonlinearity support for CW. In the future, I decide to make the index z-dependent, which is commonly-known as "wave propagation method". However, compared to them, this model will have a full nonlinear support in both solid and gas, as well as supporting pulsed scenarios.
* 3/15/2025:<br>
Finally finished implementing photoionization in both gases and solids. Current supported solid is silica only.