function PPM=resect(m,M)
%RESECT Camera resection (calibration) from 2D-3D point correspondences
%   PPM=resect(m,w) returns the camera matrix PPM given a list of (at least 6)
%   3D points w and the corresponding 2D image points m.
%
%   See also: PROJ

% Algorithm: DLT 

% Author: Andrea Fusiello, 2005

%Controllo del formato dei parametri di input
[rm,cm]=size(m);
if (rm ~= 2)
  error('Le coordinate immagine devono essere cartesiane!!');
end

[rw,cw]=size(M);
if (rw ~= 3)
  error('Le coordinate mondo devono essere cartesiane!!');
end

if ((cm < 6) | (cw < 6))
  error('Per la calibrazione servono almeno 6 punti!!!');
end

% li salvo per determinare il segno di P alla fine
m1 = m(:,1);
M1 = M(:,1);

% preconditioning    
[Tm,m] = precond(m);
[TM,M] = precond(M);

u=m(1,:);
v=m(2,:);

A=[];
for i=1:size(v,2)
  A=[A;
    0, 0, 0, 0, -M(:,i)', -1, v(i)*M(:,i)', v(i);
    M(:,i)', 1, 0, 0, 0, 0, -u(i)*M(:,i)', -u(i) ];
end


% solution vector corresponding to the 
% least singular value of A
p = ns(A);

% recover the matrix
P = reshape(p,4,3)';

% postconditioning
P = inv(Tm)*P*TM;

% bisogna anche sistemare il segno di P, in modo che i punti siano
% davanti

s =  sign(sum((P * [M1;1])./ [m1;1])/3);

% Normalizzo la matrice PPM in modo che ||q3||=1
q3=P(3,1:3);
PPM=P./s*norm(q3);



%% raffinamento non-lin
%      
%  [Kc,Rc,tc] = art(P1c);
%  [p(1),p(2:4)] = quat(Rc);
%  p(5:7) = tc;
%  k = vech(Kc');
%  p(8:12) = k(1:5);
%  
%  _resect_cost(p,m1,M);
%  
%  p = fmins('_resect_cost',p,[],[],m1,M);
%  
%  Kc=[p(8), p(9), p(10)
%    0,    p(11), p(12)
%    0,     0,       1];
%  P1c = Kc * [iquat(p(1),p(2:4)), p(5:7)'];

% 
% function c = resect_cost(p,m,w)
% 
% K = ivech([p(8:12);1])';
% 
% P = K * [iquat(p(1),p(2:4)), p(5:7)'];
% 
% mp = proj(P,w);
% 
% c = norm(m-mp,'fro')^2;
% 


