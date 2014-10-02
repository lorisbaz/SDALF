function [vo1,vo2,vo3] = absolute(X,Y,method);
%ABSOLUTE Solve absolute orientation
%
%   G = absolute(X,Y) computes the rigid transformation between two 
%   3D point sets X and Y. When G is applied to Y, it brings Y on X. 
%   Entries containing NaN in X are discarded together with the 
%   corresponding entries in Y.
%   G = absolute(X,Y,'scale') compute absolute orientation with scale.
%   
%   See also RIGID.

% Algorithm ref.: Kanatani
%
% Author: A. Fusiello, 1998

X=X';
Y=Y';

% discard NaN entries in X and correspondingly in Y
i = find(~isnan(X));
X = reshape(X(i),length(X(i))/3,3);
Y = reshape(Y(i),length(Y(i))/3,3);

dime = size(Y,1);

%% compute centroids
cm = sum(Y,1)./dime;
cd = sum(X,1)./dime;

%% subtract centroids
Yb = rigid([eye(3); -cm]',Y')';
Xb = rigid([eye(3); -cd]',X')';

if nargin == 2
    method = 'noscale';
end

if strcmp(method,'scale')
  % compute scale
  a = sqrt(sum(Xb.^2,2));
  b = sqrt(sum(Yb.^2,2));
  s = (a'*b)/(a'*a);
  
elseif strcmp(method,'noscale')
    s = 1;
else
    error('metodo inesistente');
end

% apply scale
Xb = s * Xb;
cd = s * cd;


%% compute rotation
K =  Xb' * Yb;
[U,D,V]=svd(K);
S = diag([1,1,det(U*V')]);
R = U*S*V';

%%  compute traslation
t = cd' - R * cm';

%% rigid transformation
G =  [R, t];

res = rmse(X',rigid(G,Y')./s);

%solo per debug
%fprintf('residuo di absolute orientation: %0.5g \n',res);

if nargout == 1
    vo1 = G;
end

if nargout == 2
    vo1 = G;
    vo2 = s;
  end
  
if nargout == 3
    vo1 = G;
    vo2 = s;
    vo3 = res;
end
