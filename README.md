# aSTIplus: Regularized Susceptibility Tensor Imaging With Asymmetry Constraint in the Human Brain in vivo
ASTI+ is a new proposed STI reconstruction method for human brain in vivo. 

ASTI+ improves in vivo STI reconstruction by relaxing the symmetry constraint, adding isotropic susceptibility regularization to CSF region and morphology constraint to white matter region.

This repository contains implementations of Fc-aSTI model for simulated data and in vivo human brain data.

# File Descriptions
1. simulation_data: this folder constains the ground truth of simulation experiments.
2. utils: this folder contains the source code to reconstruct aSTI+ model. including:
(1)

# Preparations before aSTI+ reconstruction
1. Extract the brain mask, white matter mask and CSF mask from magnitude images using FSL Bet and FAST;
2. Preprocess phase including Laplacian-based phase unwrapping, V_SHARP background removal using STI_Suite V3.0 toolbox.
3. Coregister the magnitude images acquired at the first echo of different orientations to the magnitude image acquired at normal supine head position using FSL FLIRT, and obtain the vector of measured field shifts and the unit vectors of applied magnetic field strength in the subject frame of reference.
4. Reconstruct QSM image using STAR-QSM using STI_Suite V3.0 toolbox and derive the smoothed QSM through Gaussian smooth. The discrete gradient was obtained from the smoothed QSM image and the morphology mask was derived by thresholding.

Once the above steps are completed, you can run STIrecon.m function

# Simulation experiments
