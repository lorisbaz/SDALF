function [vo1, vo2] = proj(P,w,t)
%PROJ  calcola una trasformazione prospettica (dalle coordinate
%      3D alle coordinate pixel).
%
%m = proj(P,w) calcola una trasformazione prospettica P delle
%cordinate 3D w. E restituisce le coordinate immagine m con 
%precisione floating point.
%
%m = proj(P,w,'f') come la precedente.
%
%m = proj(P,w,'i') calcola una trasformazione prospettica P delle
%cordinate 3D w. E restituisce le coordinate immagine m intere. 

% Author: Andrea Fusiello

%Controllo del formato dei parametri di input
 if (ismpp(P) == 0)
     warning('Matrice di proiezione prospettica non adeguata!!');
 end

[rw,cw]=size(w);
if (rw ~= 3)
    error('Le coordinate mondo devono essere cartesiane!!');
end


if nargin < 3
    t = 'f';
end

wo =[w; ones(1,size(w,2))];
mo = P*wo ;
mo = mo(1:2,:)./ [mo(3,:)' mo(3,:)']'; 

switch t
    case 'f' 
        m = mo(1:2,:);
    case 'i'
        m = [round(mo(1,:))
             round(mo(2,:))];
    case 'n'
        m = mo(1:2,:) + 0.1*rand(size(mo(1:2,:)));
    otherwise
        error('Tipo non conosciuto');
end
if nargout == 1
    vo1=m;
end

if nargout == 2
    vo1 = m(1,:);
    vo2 = m(2,:);
end


