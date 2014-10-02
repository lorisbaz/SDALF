function w = intersect(PPM,m,method)
%INTERSECT  Reconstructs 3D structure by triangulation
%
%w = intersect(PPM,m) calcola la ricostruzione 3D dato in insieme di
%matrici di proiezione prospettica (PPM) e un insieme di punti
%immagine corrispondenti.

% Author: Andrea Fusiello


siz = size(m);

dim =  siz(1);
numP = siz(2);
numV = siz(3);


if (numV<2)
    error('Sono necessarie almeno 2 viste!!');
end

for i = 1:numV
    if (ismpp(PPM(:,:,i)) == 0)
        warning('Matrice di proiezione prospettica non adeguata!!');
    end
end


if (dim ~= 2)
  error('Le coordinate immagine devono essere cartesiane!!');
end

if nargin == 2
    method = 'le';
end

if strcmp(method,'lls')
    w = intersect_lin(PPM,m);
elseif strcmp(method,'midpoint')
    w = intersect_mp(PPM,m);
elseif strcmp(method,'le')
    w = intersect_le(PPM,m);
else
    error('metodo inesistente');
end

