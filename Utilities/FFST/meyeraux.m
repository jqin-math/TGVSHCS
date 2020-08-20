function y = meyeraux(x)
% meyer wavelet auxiliary function:
% v(x) = 35*x^4 - 84*x^5 + 70*x^6 - 20*x^7.
%
% INPUT:
%  x                (vector) grid points
%
% OUTPUT:
%  y				(vector) values at given points x
%
%--------------------------------------------------------------------------
% 2012-01-20, v1.0, (c) Sören Häuser

% Auxiliary function values.
p = [-20 70 -84 35 0 0 0 0];
y = polyval(p,x);

%y = 0 for x<0, y=1 for x>1
int1 = y.*(x>=0).*(x<=1);
y = int1 + (x>1);
