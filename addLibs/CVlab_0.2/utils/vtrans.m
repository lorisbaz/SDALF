function b = vtrans(a, d)
% vtrans(a,d)     vec-transpose operator
% change from a d*s by n matrix to a d*n by s matrix.

% Author: Thomas P Minka (tpminka@media.mit.edu)

s = size(a);
n = s(2);
s = s(1)/d;

b = zeros(d*n, s);

if n < s
  for i = 1:n
    b((i*d-d+1):(i*d), :) = reshape(a(:, i), d, s);
  end
else
  for i = 1:s
    b(:, i) = reshape(a((i*d-d+1):(i*d), :), d*n, 1);
  end
end
