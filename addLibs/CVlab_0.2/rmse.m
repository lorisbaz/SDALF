function e = rmse(x,y)
%RMSE    Standard deviation.
%   For vectors, RMSE(X,Y) returns the root mean square error of X with
%   respect to the true value Y.
%
%   Matrices are treated as vectors.
%   
%   RMSE(X,Y) normalizes by (N-1) where N is the sequence length.  This
%   makes RMSE(X,Y)^2 the best unbiased estimate of the variance if X
%   is a sample from a normal distribution centered in Y.
%

% note:  typical residual in many problems is  norm(x-y,'fro')^2

% Algorithm ref.: none
%
% Author: A. Fusiello, 2006


% The Frobenius norm is the square root of the sum of squares of the
% entries. The typical residual in many problems is norm(x1-x2,'fro')^2

e = norm(x-y,'fro')/(sqrt(size(x,1) * size(x,2) -1));

% for vectors this is equivalent to:
% e = sqrt(sum((x - y).^2)/(length(x)-1));
