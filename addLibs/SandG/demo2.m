% this is demo2.m
pp = [200 50]/300;
mu1 = [0 0 0];
mu2 = [-6 3 6];
mu3 = [6 6 4];
mu = [mu1' mu2' mu3'];
covar(:,:,1) = diag([9 4 1]);
covar(:,:,2) = [4    -3.2 -0.2;...
               -3.2   4   0;...
               -0.2   0   1];
covar(:,:,3) = [4    3.2  2.8;...
               3.2    4   2.4; ...
               2.8   2.4   2];
npoints = 300;
y = genmix(npoints,mu,covar,pp);
clear npoints mu mu1 mu2 mu3 covar
[bestk,bestpp,bestmu,bestcov,dl,countf] = mixtures4(y,1,25,0,1e-4,0)

