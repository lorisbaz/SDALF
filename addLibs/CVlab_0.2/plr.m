function [R,P] = plr(A);
%PLR Compute polar decomposition
%
% R = plr(A) return the closet rotation matrix to A in Frobenius norm

% Author: Andrea Fusiello

[U,D,V] = svd(A);

R = U * diag([1,1,det(U*V')]) * V';

if nargout == 2
  P = V*D*V';
end

assert(isrot(R));
