function H = hm(M,m)
%HM Computes homography from point correspondences (2D or 3D).
%
%
%H = hm(M,m) date le corrispondenze tra due piani la 
%funzione calcola l'omografia H, che li lega.
%        m = H*M

% Author: Andrea Fusiello

[rM,cM]=size(M);

if  (rM == 2)
  % omografia del piano
  H = h2m(M,m);
  elseif (rM == 3)
% omografia dello spazio
 H = h3m(M,m);
else
  error('Trasformazione non implementata');
end




















