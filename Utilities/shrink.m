function y = shrink(x,mu)
% Soft-thresholding function
% y_i = sign(x_i)*max(abs(x_i)-mu,0);
y = sign(x).*max(abs(x)-mu,0);
% y = zeros(size(x));
% y(x>mu) = x(x>mu)-mu;
% y(x<mu) = x(x<mu)+mu;