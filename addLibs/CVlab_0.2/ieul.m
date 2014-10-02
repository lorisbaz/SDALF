function [a,b]  = ieul(R,t)
%IEUL compute Euler angles (yaw,pitch,roll) given a rotation matrix 
%
%IEUL funzione che calcola gli angoli [psi,theta,phi] di Eulero data la 
%     matrice di rotazione.
%
%[a,b] = ieul(R) ritorna il vettore colonna con i 3 angoli di Eulero 
%RPY dato la matrice di rotazione R rispetto alla terna fissa.
%Nel vettore a abbiamo la soluzione per theta appartenete all'intervallo
%(-pi/2,pi/2), mentre nel vettore b abbiamo la soluzione per theta 
%appartenente all'intervallo (pi/2,3pi/2).
%
%[a,b] = ieul(R,'RPY') come la precedente.
%
%[a,b] = ieul(R,'ZYZ') ritorna il vettore colonna con i 3 angoli di Eulero 
%ZYZ dato la matrice di rotazione R rispetto alla terna corrente.
%Nel vettore a abbiamo la soluzione per theta appartenete all'intervallo
%(0,pi), mentre nel vettore b abbiamo la soluzione per theta 
%appartenente all'intervallo (-pi,0).

% Authors: A. Fusiello, M. Mancini.


if (nargin < 2)
    t='RPY';
end

%Controllo che sia una matrice di rotazione
if (~isrot(R))
    error('La matrice R non è una matrice di rotazione!!')
end

switch t
    case 'RPY'
        phi = atan2(R(2,1),R(1,1)); 
        theta = atan2(-R(3,1),sqrt(R(3,2)^2+R(3,3)^2));
        psi =  atan2(R(3,2),R(3,3));

        phi1 = atan2(-R(2,1),-R(1,1)); 
        theta1 = atan2(-R(3,1),-sqrt(R(3,2)^2+R(3,3)^2));
        psi1 =  atan2(-R(3,2),-R(3,3));
    case 'ZYZ'
        phi = atan2(R(2,3),R(1,3)); 
        theta = atan2(sqrt(R(1,3)^2+R(2,3)^2),R(3,3));
        psi =  atan2(R(3,2),-R(3,1));
        
        phi1 = atan2(-R(2,3),-R(1,3)); 
        theta1 = atan2(-sqrt(R(1,3)^2+R(2,3)^2),R(3,3));
        psi1 =  atan2(-R(3,2),R(3,1));
    otherwise
        error('Tipo di angoli non supportato');
end

        
a(3)=phi;
a(2)=theta;
a(1)=psi;
a=a';

b(3)=phi1;
b(2)=theta1;
b(1)=psi1;
b=b';

