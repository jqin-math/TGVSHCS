function [ST,Psi] = shearletTransformSpect(A,numOfScales,shearlet_spect,shearlet_arg)
% the shearlet-transform of an image A for given number of scales with
% (optional) shearlet named shearlet_spect (with optional shearlet_arg).
% Default is 'meyerShearletSpect'. 
% The return value ST contains the shearlet coefficients in a 3-d-Matrix
% where the third index indicates the respective shear (1 for lowpass,
% 2:end for the different shears and scales). The images are ordered
% ascending with the scale and within each scale counter-clockwise with the
% direction of the shear. The output Psi contains the respective shearlets
% (in the fourier domain).
%
% INPUT:
%  A                (matrix) image (or data) to transform
%  numOfScales		(int) number of scales OR
%	 				(3-d-matrix) precomputed Psi
%  shearlet_spect	(string) name of shearlet(spectrum) (optional)
%  shearlet_arg		(var) optional argument for shearlet (optional)
%
% OUTPUT:
%  ST				(3-d-matrix) shearlet transform
%  Psi				(3-d-matrix) spectrum of shearlets
%
%--------------------------------------------------------------------------
% 2012-01-20, v1.0, (c) Sören Häuser

	%% initialization
	%persistent
	persistent L;
	persistent SHEARLET_SPECT;
	persistent SHEARLET_ARG;
	persistent PSI;
	persistent NUMOFSCALES;

	%images size
	l = size(A);

	%% error trapping
	%setting default values
	if(nargin<4)
		shearlet_arg = [];
	end

	if(nargin<3)
		shearlet_spect = 'meyerShearletSpect';
	end

	%compute number of scales
	if(nargin<2)
		numOfScales = floor(0.5 * log2(max(l)));
		if(numOfScales<1)
			error('image too small!');
		end
	end

	if(nargin<1)
		error('No image given!');
	end

	%check image size
	if( min(l) <2 )
		error('A has to be an image');
	end

	if(l(1) ~= l(2))
		error('A has to be a square image')
	end

	%% precompute Psi
	% if not cached compute and cache
	if((~isscalar(numOfScales)) && ndims(numOfScales == 3))
		Psi = numOfScales;
	%all parameters the same as last time? - changed parameters need recomputation of Psi
	elseif(isempty(L) || any(l~=L) || numOfScales ~= NUMOFSCALES || ~strcmp ( shearlet_spect, SHEARLET_SPECT ) || any ( shearlet_arg ~= SHEARLET_ARG ) )
		%check for interruption
		if( L == 0 )
			warning('Either recovering from previous interruption, or current image is garbage!');
		end

		%Provoke recomputation upon interruption
		L = 0;

		%Psi
		Psi = scalesShearsAndSpectra ( l, numOfScales, shearlet_spect, shearlet_arg );

		%save persistent variables
		L = l;
		SHEARLET_SPECT = shearlet_spect;
		SHEARLET_ARG = shearlet_arg;
		PSI = Psi;
		NUMOFSCALES = numOfScales;
	else
		Psi = PSI;
	end

	%% shearlet transform
	uST = Psi .* repmat(fftshift(fft2(A)),[1,1,size(Psi,3)]);
	ST = ifft2(fftshift(fftshift(uST,1),2));
	clear uST;
		
end
