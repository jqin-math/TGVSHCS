%% TGV2-Shearlet L1 Based Image Reconstruction with ADMM demo
% Jing, 10-27-2012

%% Load the image data for test
clear;clc;close all;
path(path,genpath(pwd));
load grayimages barbara256

%% Simulate the partial Fourier data
I     = barbara256;
I     = imnormalize(I);
[m,n] = size(I);

% If the image is NOT square, squarize it.
if m~=n
    I = squarize(I);
end
N     = m*n;

b     = fft2(I);
% Radial sampling scheme
ls    = 40;
pick  = fftshift(MRImask(n,ls));
p     = double(pick);
pr    = sum(p(:))/N;
bv    = p.*b;

%% TGV2SHCS ADMM
n = size(I,1);
L = 2; % Number of shearlet scales
para.kernels       = shearletkernels_FFST(n,L);
para.trans         = @(x)(shearlet2_FFST(x,L));
para.lambda        = 1e-2;
para.beta          = 1e3;
para.alpha         = [.0008; .001];
para.mu            = [1e3; 1e-3; 1e-5];
para.gamma         = [1e-2; 1e-2; 1e-2];
para.maxiter       = 500;
para.correctopt    = 'abs';

output = TGV2L1L2_ADMM(I,p,bv,para);

%% Display the results
% Display the reconstructed image
figure
subplot(131)
imagesc(I); axis image off;colormap gray;
title('ground truth')
subplot(132)
imagesc(fftshift(p));axis image off;colormap gray;
title(sprintf('sampling rate %4.2f%%',pr*100))
subplot(133)
imagesc(output.u); axis image off;colormap gray;
title(sprintf('recon. image, rel err %4.2f%%',output.enderr*100))

% Plot the relative errors v.s. iterations
figure
plot(output.err)
title('convergence')
xlabel('iteration number')
ylabel('relative error')