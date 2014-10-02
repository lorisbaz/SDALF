function k = commutation(n, m)
% commutation(n, m) or commutation([n m])
% returns Magnus and Neudecker's commutation matrix of dimensions n by m

% Author: Thomas P Minka (tpminka@media.mit.edu)

if nargin < 2
  m = n(2);
  n = n(1);
end

if 0
  % first method
  i = 1:(n*m);
  a = reshape(i, n, m);
  j = vec(a');
  k = zeros(n*m,n*m);
  for r = i
    k(r, j(r)) = 1;
  end
else
  % second method
  k = reshape(kron(vec(eye(n)), eye(m)), n*m, n*m);
end

% vec(kron(vec(eye(n)), eye(m)))
% vtrans(kron(kron(eye(n), eye), eye(m))*kron(vec(eye(n)),vec(eye(m))), n*m)
