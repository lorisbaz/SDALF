function y = multinorm(x,m,covar)
% Evaluates a multidimensional Gaussian
% of mean m and covariance matrix covar
% at the array of points x
%
[dim npoints] = size(x);
dd = det(covar+ realmin*eye(dim));
in = inv(covar+ realmin*eye(dim));
ff = ((2*pi)^(-dim/2))*((dd)^(-0.5));
quadform = zeros(1,npoints);
centered = (x-m*ones(1,npoints));
if dim ~= 1
   y = ff * exp(-0.5*sum(centered.*(in*centered)));
else
   y = ff * exp(-0.5*in*centered.^2 );
end



