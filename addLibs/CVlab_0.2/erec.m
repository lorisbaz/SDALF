function PPM = erec(imagePoints,A)
%EREC Computes a Euclidean reconstruction (motion only) 
%  
% Algorithm: factorization of essential matrices
%
% Author: Andrea Fusiello, 2005

siz = size(imagePoints);

numberOfPoints = siz(2);
numberOfViews = siz(3);
% 


for i=1:numberOfViews
  for j=i+1:numberOfViews
    % calcola la matrice fondamentale ij 
    Fs(:,:,i,j) = fm(imagePoints(:,:,i),imagePoints(:,:,j));  
  end
  
end


% first view
PPM(:,:,1) = A *[ 
  1     0     0     0
  0     1     0     0
  0     0     1     0 ];


% second view
[R,t] = sr(Fs(:,:,1,2),A,A,imagePoints(:,:,1),imagePoints(:,:,2));

% t12 comes from the groung truth; it fixes the scale factor
% t = t/norm(t)*norm(t12);

PPM(:,:,2) = A * [R   t];
t12 = t;

% any view > 3
for i=3:numberOfViews
  
  % fattorizzazione RS  
  [R,t]     = sr(Fs(:,:,1,i),A,A,imagePoints(:,:,1),imagePoints(:,:,i));  
  [R2i,t2i] = sr(Fs(:,:,2,i),A,A,imagePoints(:,:,2),imagePoints(:,:,i));
  
  % normalization factor to obtain coherent projection matrices
  % Explanation: see Luong & Faugeras IJCV 33(3), 1997 pg 32
  % Algorithm: see Zeller & Faugeras, INRIA RR 2793, sec. 6
  k = (norm(cross(R2i*t12,t2i))^2)/(cross(t,t2i)' * cross(R2i*t12,t2i));
  
  % camera matrix 
  PPM(:,:,i) = A * [R, k*t];
  
end

% at this point, the reconstructed camera matrices are in PPM     

