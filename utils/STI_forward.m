function deltaB = STI_forward(chitensor,H0subArray)
OriNum = size(H0subArray,1); %H - N*3
sizeVol=[size(chitensor,1),size(chitensor,2),size(chitensor,3)];
deltaB=zeros([sizeVol,OriNum]);

chi_k = zeros(size(chitensor));
for i=1:6
    chi_k(:,:,:,i)=fftn(chitensor(:,:,:,i));
end

Nx = sizeVol(1);
Ny = sizeVol(2);
Nz = sizeVol(3);
[KY_Grid, KX_Grid, KZ_Grid] = meshgrid(-Ny/2:Ny/2-1,-Nx/2:Nx/2-1,-Nz/2:Nz/2-1);      % k-space grid
KX_Grid = fftshift(KX_Grid);      KY_Grid = fftshift(KY_Grid);      KZ_Grid = fftshift(KZ_Grid);
KSq = KX_Grid.^2 + KY_Grid.^2 + KZ_Grid.^2;          % k^2
clear Nx Ny Nz

a11_array = zeros([sizeVol, OriNum]);
a12_array = zeros([sizeVol, OriNum]);
a13_array = zeros([sizeVol, OriNum]);
a22_array = zeros([sizeVol, OriNum]);
a23_array = zeros([sizeVol, OriNum]);
a33_array = zeros([sizeVol, OriNum]);

for OriInd = 1:OriNum
    H = H0subArray(OriInd,:)./norm(H0subArray(OriInd,:));
    h1 = H(1); 
    h2 = H(2);
    h3 = H(3);
    
    KHdKSq = (KX_Grid*h1 + KY_Grid*h2 + KZ_Grid*h3)./KSq; % nan at the center
    KHdKSq(isnan(KHdKSq)) = 0;
    a11_array(:,:,:,OriInd) = h1.^2/3 - KHdKSq.*KX_Grid*h1;
    a22_array(:,:,:,OriInd) = h2.^2/3 - KHdKSq.*KY_Grid*h2;
    a33_array(:,:,:,OriInd) = h3.^2/3 - KHdKSq.*KZ_Grid*h3;

    a12_array(:,:,:,OriInd) = 2*h1.*h2/3 - KHdKSq.*(KX_Grid*h2 + KY_Grid*h1);
    a13_array(:,:,:,OriInd) = 2*h1.*h3/3 - KHdKSq.*(KX_Grid*h3 + KZ_Grid*h1);
    a23_array(:,:,:,OriInd) = 2*h2.*h3/3 - KHdKSq.*(KY_Grid*h3 + KZ_Grid*h2); 
end
% remove nan
a11_array(isnan(a11_array)) = 0;  
a12_array(isnan(a12_array)) = 0;   
a13_array(isnan(a13_array)) = 0;
a22_array(isnan(a22_array)) = 0;  
a23_array(isnan(a23_array)) = 0;   
a33_array(isnan(a33_array)) = 0;

clear KHdKSq KSq KX_Grid KY_Grid KZ_Grid 

disp('Simulating with given tensor ...')

for ori = 1:OriNum      
    delta_k = a11_array(:,:,:,ori).*chi_k(:,:,:,1) + a12_array(:,:,:,ori).*chi_k(:,:,:,2) + a13_array(:,:,:,ori).*chi_k(:,:,:,3) + ...
                a22_array(:,:,:,ori).*chi_k(:,:,:,4) + a23_array(:,:,:,ori).*chi_k(:,:,:,5) + a33_array(:,:,:,ori).*chi_k(:,:,:,6) ;
    deltaB(:,:,:,ori)  = real(ifftn(delta_k));    
end
end