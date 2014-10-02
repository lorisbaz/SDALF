function [vo1,vo2,vo3] = p3t(T,w)
%_P3T Applica una trasformazione proiettiva nel 3D.
%
%wt = p3t(T,w) data una trasformazione T (omogenea)
%e un insime di punti mondo w, permette di calcolare
%una trasformazione dei punti nello spazio 3D.

% Author: Andrea Fusiello

[na,ma]=size(T);
if na~=4 | ma~=4
    error('Formato errato della matrice di trasformazione (4x4)!!');
end

[rw,cw]=size(w);
if (rw ~= 3)
    error('Le coordinate mondo devono essere cartesiane!!');
end


tmp = T * [w; ones(1,size(w,2))]; 
tmp = tmp(1:3,:)./ [tmp(4,:)' tmp(4,:)' tmp(4,:)']';
wt= tmp(1:3,:);


if nargout == 1
    vo1=wt;
end

if nargout == 3
    vo1 = wt(1,:);
    vo2 = wt(2,:);
    vo3 = wt(3,:);
end
