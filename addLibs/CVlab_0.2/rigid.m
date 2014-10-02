function wr  = rigid(G,w)
%RIGID applica una trasformazione rigida.
%
%wr  = rigid(G,w) applica una trasformazione rigida G
%all'insieme di punti w. La trasformazione deve essere
%una matrice 3x4

% Author: Andrea Fusiello

[na,ma]=size(G);

if (na~=3 & na~=4) | ma~=4
    error('Formato errato della matrice di trasformazione!!');
end


if na==4| ma==4
    G = G(1:3,:);
end


[rw,cw]=size(w);
if (rw ~= 3)
    error('Le coordinate mondo devono essere cartesiane!!');
end


HM =[w; ones(1,size(w,2))];
wr = (G*HM);
