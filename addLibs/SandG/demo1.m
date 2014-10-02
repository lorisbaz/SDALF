% this is demo1.m
pp = [0.3333 0.3333];
mu1 = [0 0];
mu2 = [0 2];
mu3 = [0 -2];
mu = [mu1' mu2' mu3' ];
covar(:,:,1) = [2 0; 0 0.2];
covar(:,:,2) = [2 0; 0 0.2];
covar(:,:,3) = [2 0; 0 0.2];
y = genmix(900,mu,covar,pp);
clear covar mu mu1 mu2 mu3
[bestk,bestpp,bestmu,bestcov,dl,countf] = mixtures4(y,1,25,0,1e-4,0)
