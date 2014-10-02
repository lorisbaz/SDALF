function H = h2m(M,m)
%_H2M Calcola l'omografia tra due immagini
%
%H = h2m(M,m) date le corrispondenze tra due piani la 
%funzione calcola l'omografia H, che li lega.
%        m = H*M
%
%N.B.:
%
%Sono necessari almeno 4 punti.

%    Algorithm: DLT (see Zisserman for example)
%
%    Author: A. Fusiello 2000


[rM,cM]=size(M);
if (rM ~= 2)
    error('Le coordinate immagine devono essere cartesiane!!');
end

[rm,cm]=size(m);
if (rm ~= 2)
    error('Le coordinate immagine devono essere cartesiane!!');
end


n_points = size(M(1,:),2);
if n_points < 4
  error('Sono necessari almeno 4 punti!!!');
end


%pre-conditioning 
 
[TM,M] = precond(M);
[Tm,m] = precond(m);


% $$$ % scale the point coordinates so that centroid is in origin 
% $$$ % and average radius  is sqr(2)
% $$$ xc = (sum(M(1,:))/n_points);
% $$$ yc = (sum(M(2,:))/n_points);
% $$$ uc = (sum(m(1,:))/n_points);
% $$$ vc = (sum(m(2,:))/n_points);
% $$$ M(1,:) = M(1,:) - xc;
% $$$ M(2,:) = M(2,:) - yc;
% $$$ m(1,:) = m(1,:) - uc;
% $$$ m(2,:) = m(2,:) - vc;
% $$$ rM = sum(sqrt(M(1,:).^2 + M(2,:).^2))/n_points/sqrt(2);
% $$$ rm = sum(sqrt(m(1,:).^2 + m(2,:).^2))/n_points/sqrt(2);
% $$$ M(1,:) = M(1,:)/rM;
% $$$ M(2,:) = M(2,:)/rM;
% $$$ m(1,:) = m(1,:)/rm;
% $$$ m(2,:) = m(2,:)/rm;
% $$$ TM = [1 0 -xc; 
% $$$       0 1 -yc; 
% $$$       0 0  rM];
% $$$ Tm = [1 0 -uc; 
% $$$       0 1 -vc; 
% $$$       0 0  rm];
% $$$   
% $$$ % end of scaling

% Preparing the equation matrix; 
% D=[];
% z=[];
% for i=1:size(m,2)
%    D=[D;
%       M(1,i), M(2,i) , 1, 0, 0, 0, -M(1,i)*m(1,i), -M(2,i)*m(1,i);
%       0,0,0,M(1,i), M(2,i) , 1,   -M(1,i)*m(2,i), -M(2,i)*m(2,i)];
%    z = [z;
%         m(1,i); 
%         m(2,i)];
% end
% % done
% 
% % solve with pseudoinverse
% x = pinv(D)*z;
% x=[x;1];
% H = reshape(x,3,3)';
% 
% % apply the inverse scaling ...
% H = inv(Tm)*H*TM;               
% 
% % ... and normalize
% H=H./H(3,3)



A=[];
for i=1:size(m,2)
   A=[A;
      0, 0, 0, -M(1,i), -M(2,i), -1, m(2,i)*M(1,i), m(2,i)*M(2,i), m(2,i);
      M(1,i),   M(2,i),  1,  0, 0, 0,  -m(1,i)*M(1,i), -m(1,i)*M(2,i), -m(1,i)];
end
% done



% solution vector corresponding to the 
% least singular value of A
h = ns(A);

% recover the matrix
H = reshape(h,3,3)';

% apply the inverse scaling ...
H = inv(Tm)*H*TM;               

% ... and normalize
H=H./H(3,3);


