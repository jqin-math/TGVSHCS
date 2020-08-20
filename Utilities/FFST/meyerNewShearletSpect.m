function Psi = meyerNewShearletSpect( x,y,a,s,shearlet_arg )
% returns the spectrum of the shearlet "meyerShearlet" with new smooth 
% construction for given scale a, shear s, and grid xi_x and xi_y. 
% shearlet_arg is optional.
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
			Psi = meyerScaling(x) .* meyerScaling(y);
			return;
		end
	end

	%compute scaling and shearing
	asy = s .* sqrt(a) .* x + sqrt(a) .* y;
	y = a .* y;
	x = a .* x;

	%set values with x=0 to 1 (for division)
	xx = (abs(x)==0) + (abs(x)>0).*x;
	
	%compute spectrum
	W = sqrt((meyerScaling(2^(-2)*x).*meyerScaling(2^(-2)*y)).^2- (meyerScaling(x).*meyerScaling(y)).^2);
	Psi = W .* bump(asy./xx);	
end

