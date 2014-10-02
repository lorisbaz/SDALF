function Ris = ismpp(P)
%ISMPP indica se una matrice è una matrice di proiezione
%      prospettica.
%
%Ris = ismpp(P) controlla se la matrice P è una matrice
%di proiezione prospettica adeguata. In caso affermativo
%ritorna 1, altrimenti ritorna 0.

% Author: Michele Mancini

[n,m]=size(P);
if ((n ~= 3) | (m ~= 4 ))
    Ris = 0;
    return;
end


if rank(P)~=3
    Ris = 0;
    return;
end

Ris= 1;
