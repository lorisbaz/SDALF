function v = ns(A);
%NS Solve the null-space problem A*v=0 
%
% NS return the 1-d null-space of A 

% Author: Andrea Fusiello

[U, D, V] = svd(A);
D = diag(D);
% $$$ % set a tolerance *much* higher than default for rank
% $$$ tol = max(size(A))* D(1) * 1e-3;
% $$$ nullity = size(A,2) - rank(A,tol);
% $$$ % check the nullity of A (dim. ker(A)); must be 1
% $$$ if (nullity) ~= 1
% $$$   warning('ns: nullity of A is %d',nullity);
% $$$ end

% check condition number
c = D(1)/D(end-1);
if (c>200)
  warning('ns: condition number is %0.f',c);
end


v = V(:,end);


