clear
close all

path = 'Z:\yuting\STI\sim\ICBM152_norm2\'; % path of phase data
save_path = 'Z:\yuting\STI\sim\ICBM152_norm2\17vecs\'; % path to save STI maps
mkdir(save_path)

%% parameter preparation
wm = load_untouch_nii([path,'wm.nii']);
wm = wm.img;
CSF = load_untouch_nii([path,'CSF.nii']);
CSF = CSF.img;
brain_mask = load_nii([path,'mask.nii']); 
brain_mask = brain_mask.img;
morphology_mask = load_nii([path,'swG05_noise005.nii']);
morphology_mask = morphology_mask.img;
load([path,'30_phi17_noise005.mat']); 
load([path,'H17_Matrix_30.mat']);
STIParams.WMMask = wm; % white matter mask
STIParams.CSFMask = CSF; % CSF mask
STIParams.BrainMask = brain_mask; % brain mask
STIParams.PhaseImage = phase_tissue;% normalized and regularized magnetic field shift images
STIParams.H0subArray = H_Matrix; % the unit vector of the applied main magnetic field in the subject frame of reference 
STIParams.wG = morphology_mask; % morphologic mask
STIParams.sizeVol = [182,218,182]; % matrix size of magnetic field shift image
STIParams.OriNum =17; % the number of head orientations
voxel_size = [1,1,1]; % voxelsize

%% Model solution parameter setting
maxit = 1000; % maximum number of iterations 
tol = 5e-3; % LSQR tolerance (5e-3 for in vivo data and simulated data with 5% noise; 1e-3 for simulated data with 1% noise)
alpha = 3;
beta  = 1;

%% solve STI with aSTI+ model
[chi11, chi12, chi13, chi21, chi22, chi23, chi31, chi32, chi33, flag, relres, iter, resvec] = aSTIplus(STIParams, maxit, tol,alpha,beta);

save([save_path,'chi_tensor.mat'],'chi11','chi12','chi13','chi21','chi22','chi23','chi31','chi32','chi33');

%% calculate parametric descriptiopns
[MMS, MSA, cMSA, PEV, abs_PEV, chitensor, symm_part] = stimap(chi11, chi12, chi13, chi21, chi22, chi23, chi31, chi32, chi33,STIParams.sizeVol);

save_nii(make_nii(MMS,voxel_size),[save_path,'MMS.nii']);
save_nii(make_nii(MSA,voxel_size),[save_path,'MSA.nii']);
save_nii(make_nii(cMSA,voxel_size),[save_path,'cMSA.nii']);
save_nii(make_nii(abs_PEV,voxel_size),[save_path,'abs_PEV.nii']);
save([save_path,'PEV.mat'],'PEV');