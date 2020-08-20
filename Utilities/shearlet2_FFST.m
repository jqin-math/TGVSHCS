function [Y, band] = shearlet2_FFST(X,L)
% 2D shearlet transform based on FFST toolbox
% X: input 2D image matrix
% L: level of shearlet kenels
% kenels: 3D matrix obtained by scalesShearsAndSpectra(,,)
% Jing, 2-27-2012

shear = 'meyerShearletSpect';
arg   = [];
% kenels = scalesShearsAndSpectra(size(X),L,shear,[]);
% 
% K = size(kenels,3);
% Y = repmat(zeros(size(X)),[1,1,K]);
% 
% band = cell(K,1);
% 
% for i = 1:K
%     temp = fftshift(kenels(:,:,i)).*fft2(X);
%     Y(:,:,i) = ifft2(temp);
%     band{i} = Y(:,:,i);
% end

Y = shearletTransformSpect(X,L,shear,arg);
K = size(Y,3);
band = cell(K,1);
for i = 1:K
    band{i} = Y(:,:,i);
end
