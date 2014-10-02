function G = lowe(A,model3d,data2d,G0);
%LOWE solve exterior orientation with Lowe's algorithm
%
%    G = exterior(A,model3d,data2d,G0) returns camera pose G given a list of
%    3D points (model3d), the corresponding 2D image points (data2d) and the
%    intrinsic parameters (A). Being an iterative algorithms, it requires a
%    starting guess (G0).
%
%    See also: RESECT

% Algorithm ref.: Lowe

% Author: A. Fusiello, 1998, 2004

%Controllo del formato dei parametri di input
[nA,mA]=size(A);
if nA~=3 | mA~=3
    error('Formato errato della matrice dei parametri intrinseci (3x3)!!')
end

[nA,mA]=size(G0);
if nA~=3 | mA~=4
    error('Formato errato della matrice dei parametri estrinseci (3x4)!!')
end

[rml,cml2]=size(data2d);
 
if (rml ~= 2)
  error('Le coordinate 3D  devono essere cartesiane!!');
end


if (cml2<6)
  error('Ci devono essere almeno 6 punti corrispondenti!!');
end

[rml,cml3]=size(model3d);

if (rml ~= 3)
  error('Le coordinate immagine devono essere cartesiane!!');
end

if (cml3<6)
  error('Ci devono essere almeno 6 punti corrispondenti!!');
end

if (cml2 ~= cml3)
  error('i punti 2d e 3d devono essere in corrispondenza')
end


a = [ieul(G0(1:3,1:3)) ; G0(:,4)];

% Il vettore a e' composto da 6 elementi:
% i primi 3 sono gli angoli di eulero che specificano
% rotazioni attorno ad x, y, e z (in questo ordine) 
% e gli ultimi 3 corrispondono al vettore traslazione.

% constants
vanishing = 1e-15;
maxiter = 300;

% initialization
res = 1E40 * ones(2*size(data2d,2),1);
iter = 0;

% change to normalized coordinates
% (or image coordinates) 
m = pt(inv(A),data2d);

% stack point coordinates 
p0 = vec(m);

while (norm(res)>1e-12 & iter<maxiter)
    
    % proietta model3d secondo a
    p = fp(a,model3d');
    
    % update residual
    prevres=res;
    res = (p-p0);
    
    % calcola lo jacobiano del residuo
    J = jfp(a,model3d');
    
    % calcola incremento con min. quad.
    da = - pinv(J) * res;
    
    % aggiorna la trasf. rigida
    a = a+da;
    
    % condizione di uscita
    if (norm(res - prevres) < vanishing)  
        break; end
    iter = iter +1;
    
end

if (iter == maxiter) warning('raggiunto numero massimo di iterazioni');
  norm(res - prevres)
end

% risultato
G = [eul(a(1:3)) a(4:6)];


%-------------------------------------
function [p] = fp(a,c3d)
%exterior: projects 3D points according to extrinsic params in a 

P = [eul(a(1:3)) a(4:6)];

h3d =[c3d ones(size(c3d,1),1)]';
h2d = P*h3d ;
c2d = h2d(1:2,:)./ [h2d(3,:)' h2d(3,:)']'; 
u = (c2d(1,:))';
v = (c2d(2,:))'; 

%%p = matrix([u,v]',1,2*size([u,v],1))';
p = vec([u,v]');


%-------------------------------------
function J  = jfp(a,c3d)
%%exterior: compute jacobian of func

G = [eul(a(1:3)) a(4:6)];
h3d =[c3d ones(size(c3d,1),1)]';
h2d = G*h3d ;

x = h2d(1,:)';
y = h2d(2,:)';
z = h2d(3,:)';

dime = size(x,1);

J =[];
for i = 1:dime
    
    col = [ -x(i)*y(i),     -z(i)^2-y(i)^2
        z(i)^2+x(i)^2,  x(i)*y(i)    
        -y(i)*z(i),     x(i)*z(i)      
        z(i),           0   
        0	,             z(i)
        -x(i) ,         -y(i)]./(z(i)^2);
    
    J = [J col];  
end
J=J';


