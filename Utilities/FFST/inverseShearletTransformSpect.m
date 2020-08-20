function A = inverseShearletTransformSpect (ST, Psi, shearlet_spect, shearlet_arg)
% computes the inverse shearlet transform for given sheaerlet coefficients.
% If the shearlets are not given they are computed. shearlet_spect and
% shearlet_arg are optionally.
%
% INPUT:
%  ST				(3-d-matrix) shearlet transform
%  Psi				(3-d-matrix) spectrum of shearlets (optional)
%  shearlet_spect	(string) name of shearlet(spectrum) (optional)
%  shearlet_arg		(var) shearlet argument (optional)
%
% OUTPUT:
%  A				(matrix) retrieved image
%
%--------------------------------------------------------------------------
% 2012-01-20, v1.0, (c) Sören Häuser

	%% initialization / error trapping

	if(nargin<1)
		error('No coefficients given!');
	end

	if(ndims(ST) ~= 3)
		error('Wrong number of shearlet coefficients');
	end

	%setting default values
	if(nargin<3)
		shearlet_spect = 'meyerShearletSpect';
	end

	if(nargin<4)
		shearlet_arg = [];
	end

	%no shearlet given -> computation possible
	if(nargin<2)
		%compute necessary data out of ST

		%image size
		l = [size(ST,1),size(ST,2)];

		%number of scales
		%possible: 1, 4, 8, 16, 32,...
		% -> -1 for lowpass
		% -> divide by for (1, 2, 4, 8, ...
		% -> +1 results in a 2^# number -> log returns #
		numOfScales = log2((size(ST,3)-1)/4 + 1);

		%compute shearlets
		Psi = scalesShearsAndSpectra (l, numOfScales,  shearlet_spect, shearlet_arg );
	end

	%% inverse shearlet transform
	A = sum(fftshift(fftshift(fft2(ST),1),2).*Psi,3);
	A = real(ifft2(ifftshift(A)));

end
