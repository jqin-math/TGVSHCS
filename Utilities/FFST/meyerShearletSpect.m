function Psi = meyerShearletSpect( x,y,a,s,shearlet_arg )
% returns the spectrum of the shearlet "meyerShearlet" for given scale a,
% shear s, and grid xi_x and xi_y. shearlet_arg is optional.
%
% INPUT:
%  x	(meshgrid) the meshgrid for the x-axis
%  y	(meshgrid) the meshgrid for the y-axis
%  a    (real) scale
%  s    (real) shear
%  shearlet_arg	(var) optional argument for shearlet
%
% OUTPUT:
%  Psi	(matrix) spectrum
%
%--------------------------------------------------------------------------
% 2012-01-20, v1.0, (c) Sören Häuser

	if(nargin>4)
		if(strcmp(shearlet_arg,'scaling'))
			Psi = meyerScalingSpect(x);
			return;
		end
	end

	%compute scaling and shearing
	y = s .* sqrt(a) .* x + sqrt(a) .* y;
	x = a .* x;

	%set values with x=0 to 1 (for division)
	xx = (abs(x)==0) + (abs(x)>0).*x;
	
	%compute spectrum
	Psi = meyerWavelet(x).*bump(y./xx);
end
