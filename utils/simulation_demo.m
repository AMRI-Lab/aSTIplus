clear;

%% phase simulation
save_path = '../simulation_exp/'; % please define the path to save results
mkdir(save_path);

load('../simulation_data/GT_stitsr.mat');
load('../simulation_data/mask.mat');

OriNum = 17;
rot_angle = 30;
repeat_num = 5;

simulate_phase(GT_stitsr, mask, save_path, OriNum, rot_angle, repeat_num)

%% morphology mask calculation (take 5% noisy phase as an example)
n = 1;
load([save_path,num2str(n),'_phi',num2str(OriNum),'_noise005.mat']);

% QSM recon
gamma = 42.58;
B0 = 3;
TE = 0.0225; % in the unit of s
H = [0 0 1];
phase_tissue = phase_tissue(:,:,:,1)*gamma*B0*TE*2*pi;
QSM = QSM_star(phase_tissue,mask,'TE',TE*1000,'B0',B0,'H',H,'padsize',[0 0 0],'voxelsize',[1 1 1]);

% Gaussian smooth
s_QSM = imgaussfilt(QSM,2);
s_QSM = s_QSM.*mask;
percentage = 0.5; % to contain 50% of the image voxels in the brain
wG = gradient_mask_all(s_QSM, mask, percentage);
save_nii(make_nii(double(wG)),[save_path,'swG05_noise005.nii'])

%% STI recon 
load('../simulation_data/wm.mat');
load('../simulation_data/CSF.mat'); 
load('../simulation_data/mask.mat');
STIParams.WMMask = wm; % white matter mask
STIParams.CSFMask = CSF;% CSF mask
STIParams.BrainMask = mask; % brain mask

path = save_path;
load([path,num2str(n),'_phi',num2str(OriNum),'_noise005.mat']);
load([path,num2str(n),'_H',num2str(OriNum),'_Matrix_',num2str(rot_angle),'.mat']);
wG = load_nii([path,'swG05_noise005.nii']);
wG = wG.img;
STIParams.PhaseImage = phase_tissue; % normalized and regularized magnetic field shift images
STIParams.H0subArray = H_Matrix; % the unit vector of the applied main magnetic field in the subject frame of reference 
STIParams.wG = wG; % morphologic mask
STIParams.sizeVol = [182,218,182]; % matrix size of magnetic field shift image
STIParams.OriNum =17; % the number of head orientations
voxel_size = [1,1,1]; % voxelsize

% Model solution parameter setting
maxit = 1000; % maximum number of iterations 
tol = 5e-3; % LSQR tolerance (5e-3 for in vivo data and simulated data with 5% noise; 1e-3 for simulated data with 1% noise)
alpha = 3;
beta  = 1;

% solve STI with aSTI+ model
[chi11, chi12, chi13, chi21, chi22, chi23, chi31, chi32, chi33, flag, relres, iter, resvec] = aSTIplus(STIParams, maxit, tol,alpha,beta);

save([save_path,'chi_tensor.mat'],'chi11','chi12','chi13','chi21','chi22','chi23','chi31','chi32','chi33');

% parametric descriptiopns calculation 
[MMS, MSA, cMSA, PEV, abs_PEV, chitensor, symm_part] = stimap(chi11, chi12, chi13, chi21, chi22, chi23, chi31, chi32, chi33,STIParams.sizeVol);

save_nii(make_nii(MMS.*mask,voxel_size),[save_path,'MMS.nii']);
save_nii(make_nii(MSA.*mask,voxel_size),[save_path,'MSA.nii']);
save_nii(make_nii(cMSA.*mask,voxel_size),[save_path,'cMSA.nii']);
save_nii(make_nii(abs_PEV.*mask,voxel_size),[save_path,'abs_PEV.nii']);
save([save_path,'PEV.mat'],'PEV');

%% Quantitative metrics calculation
% load ground truth
load('../simulation_data/GT_MSA.mat');
load('../simulation_data/GT_MMS.mat');
load('../simulation_data/GT_PEV.mat');
wm_mask = wm;

[psnr_MSA, mssim_MSA, psnr_MMS, mssim_MMS, AE, mean_AE]=compute_metrics(GT_MSA, GT_PEV, GT_MMS, wm_mask, MSA, PEV, MMS);