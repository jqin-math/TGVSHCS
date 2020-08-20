function y = meyerWavelet(x)
% compute Meyer Wavelet
%
% INPUT:
%  x                (vector) grid points
%
% OUTPUT:
%  y				(vector) values at given points x
%
%--------------------------------------------------------------------------
% 2012-01-20, v1.0, (c) Sören Häuser

	y = sqrt(abs(meyerHelper(x)).^2+abs(meyerHelper(2*x)).^2);

end

%helper function
function y = meyerHelper(x)

	xa = abs(x);

	int1 = ((xa >= 1) & (xa < 2));
	int2 = ((xa >= 2) & (xa < 4));

	psihat = int1 .* sin(pi/2*meyeraux(xa-1));
	psihat = psihat + int2 .* cos(pi/2*meyeraux(1/2*xa-1));

	y = psihat;
end