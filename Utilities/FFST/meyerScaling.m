function y = meyerScaling(x,shearlet_arg)
% mother scaling function for meyer shearlet
%
% INPUT:
%  x                (vector) grid points
%
% OUTPUT:
%  y				(vector) values at given points x
%
%--------------------------------------------------------------------------
% 2012-01-20, v1.0, (c) Sören Häuser

	xa = abs(x);

	% Compute support of Fourier transform of phi.
	int1 = ((xa < 1/2));
	int2 = ((xa >= 1/2) & (xa < 1));

	% Compute Fourier transform of phi.
	phihat = int1 .* ones(size(x));
	phihat = phihat + int2.* cos(pi/2*meyeraux(2*xa-1));

	y = phihat;
end 