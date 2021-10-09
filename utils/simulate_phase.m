function [] = simulate_phase(GT_stitsr, mask, save_path, OriNum, rot_angle, repeat_num)
%=================================================================================
%This function is used to generate simulated data from the ground truth STI
%Inputs: 
%GT_stitsr: the ground truth susceptibility tensor (N1*N2*N3*6);
%mask: brain mask;
%save_path: path to save the generated data;
%OriNum: the number of head orientations;
%rot_angle: the largest rotation angle (in the unit of degrees);
%repeat_num: the number of repeats.
%-----------------------------------
%Outputs: 
%cMSA: MSA-weighted PEV
%chitensor: susceptibility tensor (sizeVol*9)
%symm_part: the decomposed symmetric part of the obtained asymmetric tensor(sizeVol*9)
%=====================================================================================

%% Multi Orientation Simulation
for n=1:repeat_num 
    H_Matrix=zeros(OriNum,3); 
    H_Matrix(1,:)=[0,0,1];
    %% rotation angle simulation
    for ori=2:OriNum
        alpha=rot_angle/180*pi;
        theta_x=rand(1)*2*alpha-alpha;%in-vivo totation angle between(-30°,30°)
        theta_y=rand(1)*2*alpha-alpha;
        Mrot=[1,0,0;0,cos(theta_x),-sin(theta_x);0,sin(theta_x),cos(theta_x)]*[cos(theta_y),0,sin(theta_y);0,1,0;-sin(theta_y),0,cos(theta_y)]; %rotation matrix
        temp=Mrot'*[0,0,1]';
        H_Matrix(ori,:)=temp';
    end
    save([save_path,num2str(n),'_H',num2str(OriNum),'_Matrix_',num2str(rot_angle),'.mat'],'H_Matrix');
    
    %% magnetic shift simulation 
    perturb=STI_forward(GT_stitsr,H_Matrix);
    save([save_path,'sim_deltaB',num2str(OriNum),'_',num2str(rot_angle),'.mat'],'perturb');
    
    %% add Gaussion noise
    phase_tissue=zeros(size(perturb));
    % noise level = 5%
    for i=1:size(phase_tissue,4)
        temp=perturb(:,:,:,i);
        noise_level=0.05;
        Noise=max(temp(:))*noise_level*randn(size(temp));
        phase_tissue(:,:,:,i)=perturb(:,:,:,i)+Noise;
    end
    phase_tissue=phase_tissue.*mask;
    save([save_path,num2str(n),'_phi',num2str(OriNum),'_noise005.mat'],'phase_tissue');
    
    % noise level = 5%
    for i=1:size(phase_tissue,4)
        temp=perturb(:,:,:,i);
        noise_level=0.01;
        Noise=max(temp(:))*noise_level*randn(size(temp));
        phase_tissue(:,:,:,i)=perturb(:,:,:,i)+Noise;
    end
    phase_tissue=phase_tissue.*mask;
    save([save_path,num2str(n),'_phi',num2str(OriNum),'_noise001.mat'],'phase_tissue');
end
