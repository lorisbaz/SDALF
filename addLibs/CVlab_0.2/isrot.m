function a = isrot(R)
%ISROT funzione che indica se la matrice e' una
%    matrice di rotazione.
%
%a=ISR(R) ritorna 1 se la matrice R e' una matrice
%di rotazione 0 altrimenti.

% Author: Michele Mancini

    [n,m]=size(R);
    if ((n ~= 3) | (m ~= 3 ))
        disp('La matrice di rotazione R deve essere una 3x3!!');
        a=0;
        return;
    end

    
    I = R*R';
    
    % controllo ortogonalita
    if abs(trace(I) - 3) > 10^(-12)  | abs(norm(I) - 1) > 10^(-12) 
        disp('La matrice  non è ortogonale');
        a = 0
        return
    end
    
    
    %Controllo che la norma sia unitaria
    if  abs(det(R)-1) > 10^(-12)
        disp('Il determinante della matrice deve essere pari a 1!!');
        a=0;
        return;
    end

    a=1;
