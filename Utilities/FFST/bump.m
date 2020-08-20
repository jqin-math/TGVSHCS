function y = bump(x)
% compute the function psi_2^ at given points x
%
% INPUT:
%  x                (vector) grid points
%
% OUTPUT:
%  y				(vector) values at given points x
%
%--------------------------------------------------------------------------
% 2012-01-20, v1.0, (c) Sören Häuser

y = meyerBump(1+x).*(x<=0) + meyerBump(1-x).*(x>0);
y = sqrt(y);

end

function y = meyerBump(x)

int1 = meyeraux(x).*(x>=0).*(x<=1);
y = int1 + (x>1);

end