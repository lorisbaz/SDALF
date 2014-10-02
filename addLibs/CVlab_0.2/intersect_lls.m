function w = intersect_lls(PPM,m)
%_INTERSECT_LLS Intersect with Linear Least-Squares algorithm
%
%w = intersect(PPM,m) calcola la ricostruzione 3D dato in insieme di
%matrici di proiezione prospettica (PPM) e un insieme di punti
%immagine corrispondenti.

% Algorithm: Linear LS

% Author: Andrea Fusiello

%fprintf('LIN\n');

siz = size(m);

numP = siz(2);
numV = siz(3);

w = [];

for i = 1:numP
    
    A=[];
    b=[];
    
    for view = 1:numV
        b = [b
              -PPM(1,4,view)+m(1,i,view)*PPM(3,4,view);
              -PPM(2,4,view)+m(2,i,view)*PPM(3,4,view)];
        A =  [A
              PPM(1,1:3,view)-m(1,i,view)*PPM(3,1:3,view);
              PPM(2,1:3,view)-m(2,i,view)*PPM(3,1:3,view)];
    end
    
    %cond(A)
    wp= A \ b;
    w =[w wp];
end




      
              
