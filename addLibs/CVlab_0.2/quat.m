function   [n,e] = quat(R);
%QUAT Calcola il quaternione corrispondente ad una
%     generica matrice di rotazione.
%
%[n,e] = quat(R) calcola il quaternione corrispondente 
%ad una generica matrice di rotazione restituendo la 
%parte scalare n, e la parte vettoriale e.

% Author: Andrea Fusiello

%Controllo che sia una matrice di rotazione
if (~isrot(R))
    error('La matrice R non è una matrice di rotazione!!')
end

%parte scalare del quaternione
n = 0.5*sqrt(trace(R)+1);

%parte vettotiale del quaternione
e(1) = (R(3,2)-R(2,3))/(4*n);
e(2) = (R(1,3)-R(3,1))/(4*n);
e(3) = (R(2,1)-R(1,2))/(4*n);
e=e';
