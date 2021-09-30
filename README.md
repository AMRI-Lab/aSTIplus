# aSTIplus: Regularized Susceptibility Tensor Imaging With Asymmetry Constraint in the Human Brain in vivo
ASTI+ is a new proposed STI reconstruction method for human brain in vivo. 

ASTI+ improves in vivo STI reconstruction by relaxing the symmetry constraint, adding isotropic susceptibility regularization to CSF region and morphology constraint to white matter region.

This repository contains implementations of Fc-aSTI model for simulated data and in vivo human brain data.

# Required Dependencies
1. [STI_Suite V3.0 toolbox](https://people.eecs.berkeley.edu/~chunlei.liu/software.html)
2. Nifti toolbox
3. [FSL software](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/)

# File Descriptions
1. simulation_data: this folder constains the ground truth of simulation experiments.
2. utils: this folder contains the source code to reconstruct aSTI+ model. including:  
(1) aSTIplus.m: the function of the aSTI+ model;  
(2) compute_metrics.m: the function to calculate quantitative evaluation metrics including PSNR, SSIM in MSA and MMS as well as mean angular error in white matter;  
(3) gradient_mask_all.m: the function to generate the morpology mask;  
(4) simulate_phase.m: the function to simulate phase data with different noise level from the ground truth susceptibility tensor;  
(5) simulation_demo.m: the script of the whole process of simulation experiments;  
(6) STI_forward.m: the function to generate the phase from the susceptibility tensor using STI forward model;  
(7) STI_recon.m: the script of STI reconstruction using aSTI+ method;  
(8) stimap.m: the function to generate STI parametric descriptions from the calculated asymmetric susceptibility tenssor;  

# Preparations before aSTI+ reconstruction
1. Extract the brain mask, white matter mask and CSF mask from magnitude images using FSL Bet and FAST;
2. Preprocess phase including Laplacian-based phase unwrapping, V_SHARP background removal using STI_Suite V3.0 toolbox.
3. Coregister the magnitude images acquired at the first echo of different orientations to the magnitude image acquired at normal supine head position using FSL FLIRT, and obtain the vector of measured field shifts and the unit vectors of applied magnetic field strength in the subject frame of reference.
4. Reconstruct QSM image using STAR-QSM using STI_Suite V3.0 toolbox and derive the smoothed QSM through Gaussian smooth. The discrete gradient was obtained from the smoothed QSM image and the morphology mask was derived by thresholding.

Once the above steps are completed, you can run STIrecon.m script to reconstruct susceptibility tensor

# Simulation experiments
you can run simulation_demo.m script to conduct simulation experiments, including:    
(1)generate simualted phase data with different noise level;  
(2)calculate the morphology mask from the smoothed QSM image;  
(3)reconstruct STI using aSTI+ method;  
(4)calculate quantitative metrics.
