function [psnr_MSA, mssim_MSA, psnr_MMS, mssim_MMS, AE, mean_AE]=compute_metrics(GT_MSA, GT_PEV, GT_MMS, wm_mask, MSA, PEV, MMS)
    psnr_MSA = psnr(MSA,GT_MSA);
    mssim_MSA = ssim(MSA,GT_MSA);
    
    psnr_MMS = psnr(MMS,GT_MMS);
    mssim_MMS = ssim(MMS,GT_MMS);
    
    AE = acos(abs(dot(GT_PEV,PEV,4)./(sqrt(GT_PEV(:,:,:,1).^2+GT_PEV(:,:,:,2).^2+GT_PEV(:,:,:,3).^2).*sqrt(PEV(:,:,:,1).^2+PEV(:,:,:,2).^2+PEV(:,:,:,3).^2))));
    AE(isnan(AE))=0;
    AE=AE/pi*180;
    mean_AE=mean(AE(wm_mask==1));
end