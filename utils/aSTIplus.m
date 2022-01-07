function [chi11, chi12, chi13, chi21, chi22, chi23, chi31, chi32, chi33, flag, relres, iter, resvec] = aSTIplus(STIParams, maxit, tol, alpha,beta)
OriNum = STIParams.OriNum;
sizeVol = STIParams.sizeVol;
VoxNum = prod(sizeVol);
H0subArray = STIParams.H0subArray;
deltaBArray = STIParams.PhaseImage;
BrainMask_temp = STIParams.BrainMask;
BrainMask = BrainMask_temp;
CSFMask = STIParams.CSFMask;
CSFMask = double(CSFMask);
WM_mask = STIParams.WMMask;
WM_mask = double(WM_mask);
WM_mask = WM_mask(end:-1:1,:,:);
wG = STIParams.wG;
wG = wG.*BrainMask;
TV = TVOP;
zeroMask = zeros(sizeVol(1), sizeVol(2), sizeVol(3), 'single');

%% using least square method for solving STI
% 
% setting up b for solving Ax = b;
b = zeros((OriNum+15)*VoxNum, 1, 'single');    

for OriInd = 1:OriNum    
    temp = deltaBArray(:,:,:,OriInd).*BrainMask;    
    b(((OriInd - 1)*VoxNum+1) : OriInd*VoxNum) = temp(:);
end
b((OriNum*VoxNum+1):end) = 0;           % For regularization

clear temp deltaBArray

% K space variables -----------------
Nx = sizeVol(1);
Ny = sizeVol(2);
Nz = sizeVol(3);

[KY_Grid, KX_Grid, KZ_Grid] = meshgrid(-Ny/2:Ny/2-1,-Nx/2:Nx/2-1,-Nz/2:Nz/2-1);      % k-space grid
KX_Grid = fftshift(KX_Grid);      KY_Grid = fftshift(KY_Grid);      KZ_Grid = fftshift(KZ_Grid);
KSq = KX_Grid.^2 + KY_Grid.^2 + KZ_Grid.^2;          % k^2

clear Nx Ny Nz

% precalculated coefficient matrix
a11_array = zeros([sizeVol, OriNum]);
a12_array = zeros([sizeVol, OriNum]);
a13_array = zeros([sizeVol, OriNum]);
a21_array = zeros([sizeVol, OriNum]);
a22_array = zeros([sizeVol, OriNum]);
a23_array = zeros([sizeVol, OriNum]);
a31_array = zeros([sizeVol, OriNum]);
a32_array = zeros([sizeVol, OriNum]);
a33_array = zeros([sizeVol, OriNum]);

for OriInd = 1:OriNum
    h1 = H0subArray(OriInd,1); 
    h2 = H0subArray(OriInd,2);
    h3 = H0subArray(OriInd,3);

    KHdKSq = (KX_Grid*h1 + KY_Grid*h2 + KZ_Grid*h3)./KSq; % nan at the center

    a11_array(:,:,:,OriInd) = h1.^2/3 - KHdKSq.*KX_Grid*h1;
    a22_array(:,:,:,OriInd) = h2.^2/3 - KHdKSq.*KY_Grid*h2;
    a33_array(:,:,:,OriInd) = h3.^2/3 - KHdKSq.*KZ_Grid*h3;

    a12_array(:,:,:,OriInd) = h1.*h2/3 - KHdKSq.*KX_Grid*h2;
    a13_array(:,:,:,OriInd) = h1.*h3/3 - KHdKSq.*KX_Grid*h3;
    a21_array(:,:,:,OriInd) = h2.*h1/3 - KHdKSq.*KY_Grid*h1;
    a23_array(:,:,:,OriInd) = h2.*h3/3 - KHdKSq.*KY_Grid*h3;
    a31_array(:,:,:,OriInd) = h3.*h1/3 - KHdKSq.*KZ_Grid*h1;
    a32_array(:,:,:,OriInd) = h3.*h2/3 - KHdKSq.*KZ_Grid*h2;
end

% remove nan
a11_array(isnan(a11_array)) = 0;  
a12_array(isnan(a12_array)) = 0;   
a13_array(isnan(a13_array)) = 0;
a21_array(isnan(a21_array)) = 0;  
a22_array(isnan(a22_array)) = 0;  
a23_array(isnan(a23_array)) = 0;
a31_array(isnan(a31_array)) = 0;
a32_array(isnan(a32_array)) = 0;
a33_array(isnan(a33_array)) = 0;

clear KHdKSq KSq KX_Grid KY_Grid KZ_Grid 
% ---------------------------------------------------------------

disp('solving Susceptibility Tensor Imaging (STI) with aSTI+ ...')

tic
[chi_tensor, flag, relres, iter, resvec] = lsqr(@afun,b,tol,maxit);
toc

% ------------------------------------------------------------------
% change solution format
chi11 = BrainMask.*reshape(chi_tensor(1:VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
chi12 = BrainMask.*reshape(chi_tensor((1*VoxNum+1):2*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
chi13 = BrainMask.*reshape(chi_tensor((2*VoxNum+1):3*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
chi21 = BrainMask.*reshape(chi_tensor((3*VoxNum+1):4*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
chi22 = BrainMask.*reshape(chi_tensor((4*VoxNum+1):5*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
chi23 = BrainMask.*reshape(chi_tensor((5*VoxNum+1):6*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
chi31 = BrainMask.*reshape(chi_tensor((6*VoxNum+1):7*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
chi32 = BrainMask.*reshape(chi_tensor((7*VoxNum+1):8*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
chi33 = BrainMask.*reshape(chi_tensor((8*VoxNum+1):9*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));

resvec = resvec./norm(b);                       % convert to relative residual vector

%% internal function including afun
function y = afun(x,transp_flag)
    
   if strcmp(transp_flag,'transp')              % y = A'*b
       y = zeros(9*VoxNum, 1, 'single');          % chitensor (9*Nv X 1)
       x = single(x);                             % change to single format ((Nori*Nv x 1))      
       
       for orient_i = 1:OriNum
            delta = x(((orient_i - 1)*VoxNum+1) : orient_i*VoxNum);            
            delta = reshape(delta, sizeVol(1), sizeVol(2), sizeVol(3));
            delta = fftn(delta.*BrainMask);
            
            %% from the 6N X N matrix
            chi11_temp = real(ifftn(a11_array(:,:,:,orient_i).*delta));
            chi12_temp = real(ifftn(a12_array(:,:,:,orient_i).*delta));
            chi13_temp = real(ifftn(a13_array(:,:,:,orient_i).*delta));
            chi21_temp = real(ifftn(a21_array(:,:,:,orient_i).*delta));
            chi22_temp = real(ifftn(a22_array(:,:,:,orient_i).*delta));
            chi23_temp = real(ifftn(a23_array(:,:,:,orient_i).*delta));
            chi31_temp = real(ifftn(a31_array(:,:,:,orient_i).*delta));
            chi32_temp = real(ifftn(a32_array(:,:,:,orient_i).*delta));
            chi33_temp = real(ifftn(a33_array(:,:,:,orient_i).*delta));

            y = y + [chi11_temp(:); chi12_temp(:); chi13_temp(:); chi21_temp(:);...
                chi22_temp(:); chi23_temp(:); chi31_temp(:); chi32_temp(:); chi33_temp(:)];                                    
       end
       
       % regularize off-diagonal elements in CSF to be zero
       y = y + alpha*x((OriNum*VoxNum+1):(OriNum+9)*VoxNum).*cat(1, zeroMask(:), CSFMask(:), CSFMask(:), ...
           CSFMask(:),zeroMask(:), CSFMask(:), CSFMask(:), CSFMask(:),zeroMask(:));
       
       % regularize diagonal elements in CSF to be equal
       temp = alpha*x(((OriNum+9)*VoxNum+1):(OriNum+12)*VoxNum).*repmat(CSFMask(:), [3,1]); 
       y(1:VoxNum)              = y(1:VoxNum) + temp(1:VoxNum) - temp(2*VoxNum+1:3*VoxNum);                         % x11
       y((4*VoxNum+1):5*VoxNum) = y((4*VoxNum+1):5*VoxNum) - temp(1:VoxNum) + temp(VoxNum+1:2*VoxNum);              % -x22              
       y((8*VoxNum+1):9*VoxNum) = y((8*VoxNum+1):9*VoxNum) - temp(VoxNum+1:2*VoxNum) + temp(2*VoxNum+1:3*VoxNum);   % -x33       
       
       % add morphology constraints in white matter
       temp = beta*WM_mask.*(TV'*(reshape(wG(:).*x(((OriNum+12)*VoxNum+1):(OriNum+15)*VoxNum), size(wG))));
       y(1:VoxNum) = y(1:VoxNum) + temp(:);                                 % x11
       y((4*VoxNum+1):5*VoxNum) = y((4*VoxNum+1):5*VoxNum) + temp(:);       % x22
       y((8*VoxNum+1):9*VoxNum) = y((8*VoxNum+1):9*VoxNum) + temp(:);       % x33

       disp('Iteration of transpose A ...')

   elseif strcmp(transp_flag,'notransp')            % y = A*x
        y = zeros((OriNum+15)*VoxNum, 1, 'single');        % x is chi_tensor, y is N deltaB and 15 regularization
        x = single(x);                                     % change to single format
       
        x11 = reshape(x(1:VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
        x22 = reshape(x((4*VoxNum+1):5*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
        x33 = reshape(x((8*VoxNum+1):9*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
                
        chi11k = fftn(x11);
        chi12k = fftn(reshape(x((1*VoxNum+1):2*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3)));
        chi13k = fftn(reshape(x((2*VoxNum+1):3*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3)));
        chi21k = fftn(reshape(x((3*VoxNum+1):4*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3)));
        chi22k = fftn(x22);
        chi23k = fftn(reshape(x((5*VoxNum+1):6*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3)));
        chi31k = fftn(reshape(x((6*VoxNum+1):7*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3)));
        chi32k = fftn(reshape(x((7*VoxNum+1):8*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3)));
        chi33k = fftn(x33);
       
        % data fidelity terms
       for orient_i = 1:OriNum
            delta = a11_array(:,:,:,orient_i).*chi11k + a12_array(:,:,:,orient_i).*chi12k + a13_array(:,:,:,orient_i).*chi13k + ...
                        a21_array(:,:,:,orient_i).*chi21k + a22_array(:,:,:,orient_i).*chi22k + a23_array(:,:,:,orient_i).*chi23k + ...
                        a31_array(:,:,:,orient_i).*chi31k + a32_array(:,:,:,orient_i).*chi32k + a33_array(:,:,:,orient_i).*chi33k ;

            delta = real(ifftn(delta)).*BrainMask;
           
            y(((orient_i - 1)*VoxNum+1) : orient_i*VoxNum) = delta(:);          
       end

       % regularize off-diagonal elements in CSF to be zero
       y(OriNum*VoxNum+1:(OriNum+9)*VoxNum) = alpha*x(:).*cat(1, zeroMask(:), CSFMask(:), CSFMask(:), ...
           CSFMask(:),zeroMask(:), CSFMask(:), CSFMask(:), CSFMask(:), zeroMask(:));

       % regularize diagonal elements in CSF to be equal
       y(((OriNum+9)*VoxNum+1):((OriNum+12)*VoxNum)) = alpha*repmat(CSFMask(:), [3,1]).*cat(1, ...
                        x11(:) - x22(:), x22(:) - x33(:), x33(:) - x11(:));             
             
       % add morphology constraints in white matter
       TVreg = beta*wG.*(TV*(x11 + x22 + x33)).*WM_mask;
       y(((OriNum+12)*VoxNum+1):((OriNum+15)*VoxNum)) = TVreg(:);

       disp('Iteration of nontranspose A ...')
          
   end
end

end
