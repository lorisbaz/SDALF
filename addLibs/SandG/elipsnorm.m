function y = elipsnorm(m,cov,level,dashed);
%
% draws one contour plot of a bivariate
% gaussian density with mean "m" and covariance
% matrix "cov". 
% The level is controled by "level".
% If "dashed==1", the line is dashed.
% 
if nargin<4 
   dashed=0;
end
[uu,ei,vv]=svd(cov);
a = sqrt(ei(1,1)*level*level);
b = sqrt(ei(2,2)*level*level);
theta = [0:0.01:2*pi];
xx = a*cos(theta);
yy = b*sin(theta);
cord = [xx' yy']';
cord = uu*cord;
if dashed==1
   plot(cord(1,:)+m(1),cord(2,:)+m(2),'--k','LineWidth',2)
else
   plot(cord(1,:)+m(1),cord(2,:)+m(2),'k','LineWidth',2)
end
