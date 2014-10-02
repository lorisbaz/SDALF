function S = selecth(n)
% return the matrix S such that: vech(A) = S * vec(A)
% n is the dimension of A (square)
% the matrix S is n*(n+1)/2 x n^2

S = zeros(n*(n+1)/2, n^2);

J = [1:1:n^2]';

Jm = reshape(J,n,n);

I = vech(Jm);

S(I,J) = 1;

S = diag(diag(S));


S = S(I,:);