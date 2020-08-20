function H = coneIndicator(s)
% obtain the cone indicator function for both the horizontal and the
% vertical cone
%
% INPUT:
%  s				(vector) image size
%
% OUTPUT:
%  H				(struct) H.x -> indicator for horizontal cone
%							 H.y -> indicator for vertical cone
%
%--------------------------------------------------------------------------
% 2012-01-20, v1.0, (c) Sören Häuser

leftCone = tril(ones(s)).*flipud(tril(ones(s)));

H.hor = leftCone+fliplr(leftCone);
H.ver = 1-H.hor;

end