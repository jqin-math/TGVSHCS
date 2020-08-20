function y = meyerScalingSpect(x)
% compute the spectrum of the meyer scaling function (cone adapted)
%
% INPUT:
%  x                (vector) grid points
%
% OUTPUT:
%  y				(vector) values at given points x
%
%--------------------------------------------------------------------------
% 2012-01-20, v1.0, (c) Sören Häuser

	%cone indicator
	H = coneIndicator(length(x));

	%scaling function
	yy = meyerScaling(x);

	%cone scaling
	y =  yy.*H.hor+yy'.*H.ver;
end