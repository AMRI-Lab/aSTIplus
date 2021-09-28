function [MMS, MSA, cMSA, PEV, abs_PEV, chitensor, symm_part] = stimap(chi11, chi12, chi13, chi21, chi22, chi23, chi31, chi32, chi33, sizeVol)   
%=================================================================================
%This function is used to calculate some parametric descriptions of STI
%from the obtained susceptibility tensor elements
%Inputs: 
%chi11, chi12, ..., chi33: 9 tensor elements
%sizeVol: matrix size of tensor elements
%--------------------------------------------
%Outputs: 
%cMSA: MSA-weighted PEV
%chitensor: susceptibility tensor (sizeVol*9)
%symm_part: the decomposed symmetric part of the obtained asymmetric tensor(sizeVol*9)
%=====================================================================================
    chitensor(1,1,:,:,:) = chi11;
    chitensor(1,2,:,:,:) = chi12;
    chitensor(1,3,:,:,:) = chi13;
    chitensor(2,1,:,:,:) = chi21;
    chitensor(2,2,:,:,:) = chi22;
    chitensor(2,3,:,:,:) = chi23;
    chitensor(3,1,:,:,:) = chi31;
    chitensor(3,2,:,:,:) = chi32;
    chitensor(3,3,:,:,:) = chi33;
    chitensor_T = permute(chitensor,[2,1,3,4,5]);
    
    %% tensor decomposition
    symm_part = 0.5*(chitensor_T+chitensor);
    anti_part = 0.5*(chitensor-chitensor_T);
    
    symm_part = reshape(symm_part,3,3,prod(sizeVol));
    
    %% eigenvalue and eigenvector calculation
    Chi_eig = zeros(prod(sizeVol),3);
    Chi_vec=zeros(prod(sizeVol),3,3);
    for i =1:prod(sizeVol)
        [V,D] = eig(symm_part(:,:,i));
        Chi_eig(i,:) = diag(D)';
        Chi_vec(i,:,:)=V;
    end
    Chi_eig = reshape(Chi_eig, [sizeVol, 3]);
    Chi_vec=reshape(Chi_vec,[sizeVol,3,3]);
    
    %% parametric descriptions
    MMS = mean(Chi_eig,4);
    MSA = Chi_eig(:,:,:,3) - (Chi_eig(:,:,:,1) + Chi_eig(:,:,:,2)) / 2;
    cMSA = abs(Chi_vec(:,:,:,:,3)).*repmat(MSA,[1 1 1 3]);
    abs_PEV = abs(Chi_vec(:,:,:,:,3));
    PEV = Chi_vec(:,:,:,:,3);
    
    symm_part = reshape(symm_part,[sizeVol,9]);
end