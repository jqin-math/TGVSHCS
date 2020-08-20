function H = shearletkernels_FFST(n,L,varargin)
% Generate shearlet kenels in Fourier space
shear = 'meyerShearletSpect';
arg = [];

H = scalesShearsAndSpectra([n n],L,shear,arg);
H = fftshift(fftshift(H,1),2);

if nargin>2
    k = varargin{1};
    H = H(:,:,k);
end