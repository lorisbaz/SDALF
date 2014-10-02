function y = uninorm(x,m,var)
% Evaluates a multidimensional Gaussian
% of mean m and variance var
% at the array of points x
%
ff = ((2*pi*(var+realmin))^(-1/2));
y = ff * exp((-1/(2*var))*(x-m).^2);


