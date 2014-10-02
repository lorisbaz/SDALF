function y = genmix(npoints,mean,var,prob)
%
% y = genmix(npoints,mean,var,prob)
%
% Generates 'npoints' samples from a mixture
% of Gaussian densities.
% The dimension of the space is the number of
% lines of 'mean', while the number of components, 'comp',
% of the mixture is its number of columns. That is
% each column of 'mean' is the mean of one component
% of the mixture.
% Parameter 'var' must be a 3-dimensional array such that 
% var(:,:,i) is the covariance matrix of the i-th
% component of the mixture. Of course, size(var(:,:,i))
% must be [dim,dim].
% 'probs', contains the mixing probabilities of the 
% first comp-1 components. The comp-th one is, of course,
% 1-sum(comp)
% 
%
clear y;
getpars = size(mean);
dim = getpars(1);
modes = getpars(2);
if modes~=1
   if  ( max(size(prob))~=(modes-1) |...
         min(size(prob))~=1|...
         max(prob)>1 |...
         min(prob)<0 )
      disp('Invalid vector of mixing probabilities')
      return
   end
   lastp = 1-sum(prob);
   prob = [reshape(prob,1,modes-1) lastp];
else
   prob = [1];
end
dat = []; 
for i=1:modes
   if dim~=1
      filter = sqrtm(var(:,:,i));
   else
      filter = sqrt(var(i));
   end
   len = round(npoints*prob(i));
   noi = randn(dim,len);
   dat = [dat , kron(mean(:,i),ones(1,len))+filter*noi];
end
y = dat;
