function a = invdiag_matrix(n)
% invdiag_matrix(n)
% Returns the invdiag matrix of size n.
% (labeled R_n in the paper)

% Author: Thomas P Minka (tpminka@media.mit.edu)

a = kron(eye(n), unit(1, n+1)');
a = a(1:n, 1:(n^2));
