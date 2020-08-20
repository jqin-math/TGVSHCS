function J = imnormalize(I)
% Image normalize
a = max(I(:));
b = min(I(:));
c = a-b;
J = (I-b)./c;