function y = shrink2(x,thresh)
%% Isotropic shrinkage operator in R^n spaces
% shrink2(x,thresh) returns 
%      $ y = (1-\frac{1}{thresh*\|x\|_2})x $ if $x\neq 0$
%      $ y = 0 $ otherwise.
% Inputs: 
% x - M * N * k matrix
% thresh - thresholding parameter >0.
%
% Outputs:
% y - M * N * k matrix
%
% Jing, 4-26-2013

k = size(x,3);
r = sqrt(sum(x.^2,3));
y = zeros(size(x));
for i = 1:k
    y(:,:,i) = x(:,:,i).*(1-1./r./thresh);
end
y(isnan(y)) = 0;

