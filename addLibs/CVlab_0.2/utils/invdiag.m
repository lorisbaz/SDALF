function x = invdiag(d)
% x = invdiag(d).
% inverts rdiag.

r = rows(d)/cols(d);
x = zeros(r, cols(d));
for i = 1:cols(d)
  x(:, i) = d((r*i-r+1):(r*i), i);
end
