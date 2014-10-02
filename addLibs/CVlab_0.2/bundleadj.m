function [PPMb,wb] = bundleadj(K,PPM,w,m)
% bundle adjustment with alternation, no derivatives.


siz = size(m);
numP = siz(2);
numV = siz(3);

% Precondizione: la prima tc deve essere [K|0] e la norma di w(:) unitaria

% normalize such that first camera is = [K|0] and norm(w(:))=1
G1 = [inv(K) * PPM(:,:,1); [ 0 0 0 1]];
w = p3t(G1,w);
% scale points
s = norm(w(:));
w = p3t(diag([1,1,1,s]),w);
% this is what has been applied to points:
% G1 = diag([1,1,1,s])*G1;
% apply the same (inverse) to cameras
iG1 = inv(G1)*diag([1,1,1,1/s]);
for j = 1:numV
  PPM(:,:,j) =  PPM(:,:,j) * iG1;
end 

% sanity check
assert(norm(PPM(:,:,1)- K * [eye(3) , [0 0  0]'], 'fro') < 1e-12);
assert('norm(w(:))==1');
%%fprintf('Errore di ricostruzione nelle immagini: %0.5g pixel\n',mrmse(PPM,m,w));

G = ones(size(PPM));
% obtain an initial guess; the first [RotsM|Tras] is [I|0] by default
iK = inv(K);
for i = 1:(numV-1)
  G(:,:,i) =  iK*PPM(:,:,i+1);
  % forzo la matrice di rotazione
  RotsM(:,:,i) = plr(G(1:3,1:3,i));
  Tras(:,i)  = G(1:3,4,i);
end 


%%size(RotsM)

% converts rotation matrices to quaternions
RotsQ = dcm2q(RotsM);
%% note: dcm2q returns quaternions as ROW vectors 

% put everithing in a vector (lsqnonlin wants a vector)
k=vech(K');
% Motion
Mot =  [k(1:5); vec(Tras); vec(RotsQ')];
% Stucture
Str = vec(w); 
% drop the last element because of the unit norm constraint
Str = Str(1:end-1);

iter = 0;
res1 = 0;
res2 = 10;

while(iter < 5 && abs(res1-res2)>eps)
  iter = iter + 1;
  
  % intersection
  
  % Solve with fixed motion: structure is changed  
 
  [Str,res1] = lsqnonlin(@(x) bundle_cost(x,Mot,m), Str);
   res1
  % resection
  
  % Solve with fixed structure: motion (and intrinsics) are changed
  
  [Mot,res2] = lsqnonlin(@(x) bundle_cost(Str,x,m), Mot);
  res2
 
end  

% re-build matrices from the solution vector

% intrinsic parameters
K = ivech([Mot(1:5);1])';
% translation and rotations (quaternions)
t = 3*(numV-1);
Tras = ivec(Mot(6:6+t-1),3);
RotsQ = ivec(Mot(t+6:end),4)';

% convert back form quaternions to rotation matrices 
RotsM = q2dcm(RotsQ);


% compute PPM
% the first PPM is [K|0]...
PPMb(:,:,1) = K * [eye(3), [0 0 0]'];
% ... now build the others
for i = 1:(numV-1)
  PPMb(:,:,i+1) = K * [RotsM(:,:,i), Tras(:,i)];
end

% recover last element of the structure from the unit norm constraint
Str=[Str; sqrt(1 - norm(Str)^2)];
% compute structure
wb = ivec(Str,3);


% sanity check
assert(norm(PPMb(:,:,1)- K * [eye(3) , [0 0  0]'], 'fro') < 1e-12);
assert('norm(wb(:))==1');

%%fprintf('Errore di ricostruzione nelle immagini: %0.5g pixel\n',mrmse(PPMb,m,wb));


 
%% probably it would be also necessary to precondition 3D  points?



function r = bundle_cost(Str,Mot,m)

% Str is the structure (without the first camera, which is always [K|0])
% Mot is the motion (without last element, because Y is normalized

% Author: Andrea Fusiello

siz = size(m);
numP = siz(2);
numV = siz(3);

% re-build matrices from the argument vector

% intrinsic parameters
K = ivech([Mot(1:5);1])';
% translation and rotations (quaternions)
t = 3*(numV-1);
Tras = ivec(Mot(6:6+t-1),3);
RotsQ = ivec(Mot(t+6:end),4)';

% convert back form quaternions to rotation matrices 
RotsM = q2dcm(RotsQ);


% compute PPM
% the first PPM is [K|0]...
PPM(:,:,1) = K * [eye(3), [0 0 0]'];
% ... now build the others
for i = 1:(numV-1)
  PPM(:,:,i+1) = K * [RotsM(:,:,i), Tras(:,i)];
end

% recover last element of the structure from the unit norm constraint
Str=[Str; sqrt(1 - norm(Str)^2)];
% compute structure
w = ivec(Str,3);


% compute residual vector (distances of image points)
r=[];
for i = 1:numV
  mp = proj(PPM(:,:,i),w);
  % r + r + norm(mp-m(:,:,i),'fro')^2 
  r = [r; 
       sqrt(sum((m(:,:,i)-mp).^2,1))'];

end

