function h = vech(A)
% h = vech(A)
% h is the column vector of elements on or below the main diagonal of A
% (lower triangular part of A).
% A must be square.


h = A(find(tril(ones(size(A)))));
