function  R = iquat(n,e);
%IQUAT calcola la matrice di rotazione dato il quaternione.
%
%R = iquat(n,e) calcola la matrice di rotazione dato il quaternione,
%cioe' la parte scalare n e la parte vettoriale e.

% Author: Andrea Fusiello

%Controllo che la parte vettoriale e sia un vettore colonna di 3 elementi
ne=length(e);
if (ne ~= 3)
    error('La parte vettoriale, e, deve essere un vettore di 3 elementi!!');
end


R = [ 2*(n^2+e(1)^2)-1        2*(e(1)*e(2)-n*e(3))       2*(e(1)*e(3)+n*e(2))
      2*(e(1)*e(2)+n*e(3))    2*(n^2+e(2)^2)-1           2*(e(2)*e(3)-n*e(1))
      2*(e(1)*e(3)-n*e(2))    2*(e(2)*e(3)+n*e(1))       2*(n^2+e(3)^2)-1];
