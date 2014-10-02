function PPM = calibz(m,M)
%CALIBZ Camera calibration from 2D-2D correspondences in many views
%

% Algorithm: Zhang's flexible calibration

% Author: Andrea Fusiello, 2005

siz = size(m);

numP = siz(2);
numV = siz(3);

 
    %Controllo del formato dei parametri di input
     [rm,cm]=size(m);
     if (rm ~= 2)
         error('Le coordinate immagine devono essere cartesiane!!');
     end
 
     [rw,cw]=size(M);
     if (rw ~= 2)
         error('Le coordinate modello devono essere 2D!!');
     end
 
     if ((numV < 3) )
         error('Per la calibrazione servono almeno 3  viste!!!');
     end
       
%S =  (eye(9)+commutation(3,3))*selecth(3)';
 S = duplication(3);
 

A = []; 
A1 = [];

for i =1:numV
  
  H = hm(M,m(:,:,i));
  
  HH(:,:,i) = H;
  
  A=[A
    (kron(H(:,1)',H(:,2)'))*S;
    (kron(H(:,1)',H(:,1)')-kron(H(:,2)',H(:,2)'))*S; ];
  
end


% solution vector corresponding to the 
% least singular value of A
b = ns(A);

% recover the matrix 
B = ivech(b);
B = B + B';

[vec,val]=eig(B);
% aggiusta segno di B
if sum(diag(val))<0
  val = -val;
end
% aggiusta autovalori leggermente negativi
val(val<0)=eps;
B1=vec*val*vec';

iK = real(chol(B1));
K = inv(iK);
K = K./K(3,3)

for i =1:numV
  
  H = HH(:,:,i);

  G = (1/norm(iK*H(:,1)))*iK*H;
  R = [G(:,1), G(:,2), cross(G(:,1), G(:,2))];
  t =  G(:,3);

  % forzo la matrice di rotazione 
  R = plr(R);
  
  PPM(:,:,i) =  K * [R,t];
  
end

% The condition number is large and the method is destroied even by a small
% amount of noise, even with normalization....
