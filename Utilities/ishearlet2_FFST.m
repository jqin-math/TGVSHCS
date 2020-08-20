function output = ishearlet2_FFST(X)
% 2D inverse shearlet transform based on FFST toolbox
% X: input 3D matrix 
% L: level of shearlet kernels
% kernels: 3D matrix obtained by scalesShearsAndSpectra(,,)
% Jing, 3-17-2012

shear = 'meyerShearletSpect';
if iscell(X)
    nX = size(X,1);
    Y = zeros(size(X{1},1),size(X{1},2),nX);
    for i = 1:nX
        Y(:,:,i) = X{i};
    end
    X = Y;
end

m = size(X,1); 
n = size(X,2);
L = log2((size(X,3)-1)/4+1);
K = scalesShearsAndSpectra([m,n],L,shear,[]);

output = sum(fftshift(fftshift(fft2(X),1),2).*K,3);
output = real(ifft2(ifftshift(output)));
