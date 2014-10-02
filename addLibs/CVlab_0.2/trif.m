function T  = trif(P1,P2,P3)
%TRIF calcola la matrice trifocale a partire dalle matrici 
%     di proiezione prospettica  delle tre viste.
%
% Author: Andrea Fusiello

% Controllo se sono matrici di proiezione prospettica
if ((ismpp(P1)==0) | ismpp(P2)==0 | ismpp(P3)==0)
    warning('Matrici di proiezione prospetica non adeguate!!!');
end


A2 = P2(1:3,1:3) * inv(P1(1:3,1:3));
A3 = P3(1:3,1:3) * inv(P1(1:3,1:3));

e21=P2*[-inv(P1(1:3,1:3))*P1(:,4) ;1];
e31=P3*[-inv(P1(1:3,1:3))*P1(:,4) ;1];

% calcolo matrice trifocale 
T = kron(A3,e21) - kron(e31,A2);





