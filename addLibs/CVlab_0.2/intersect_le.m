function S = intersect_le(PPM,m)
%_INTERSECT_LE Intersection with Linear-Eigen algorithm
%
%M = intersect(PPM,m) calcola la ricostruzione 3D dato in insieme di
%matrici di proiezione prospettica (PPM) e un insieme di punti
%immagine corrispondenti.

% Algorithm: Linear-Eigen

% Author: Andrea Fusiello

%fprintf('LE\n');

siz = size(m);

numP = siz(2);
numV = siz(3);


% tentativo di precondizionamento
for j = 1:numV
  [T, m(:,:,j)] = precond(m(:,:,j));
  PPM(:,:,j) = T * PPM(:,:,j);
  PPM(:,:,j) = PPM(:,:,j)./norm(PPM(3,1:3,j));
end

% auxiliary matrix containing measurements 	 
QQ = [];
for i = 1:numP
  % Q_i contiene le proiezioni dell'i-esimo punto nelle varie viste
  Q_i = [];
  for j= 1:numV
    Q_i =  [Q_i
      zeros(3,j-1) [m(:,i,j);1] zeros(3,numV-j)];
  end
  QQ(:,:,i) = Q_i;
end

% mette le matrici nella forma richiesta
P3 = [];
for j = 1:numV
  P(3*j-2:3*j,:) = PPM(:,:,j);
  P3 = [P3
    PPM(3,:,j)];
end


S = [];
for i = 1:numP  
  % triangolo il punto i
  z = ones(numV,1);
  iter =0;
  while (iter < 20)
    
    % build the weight matrix
    for j = 1:numV   
      Z(3*j-2:3*j,3*j-2:3*j) = 1/z(j)*eye(3);
    end
    
    %%fprintf('Solving point: %d\n', i );
    M = ns(Z*(QQ(:,:,i)*P3-P));
    M = M(1:3)./M(4);
    
    iter=iter+1;
    zprec = z;
    z  = P3 * [M;1];
    
    if norm(zprec-z)<1e-12
      %% fprintf('Point: %d converged\n', i );
      break;
    end
    
  end
  
  S =[S,M];
  
end


% versione vecchia, senza raffinamento iterativo

% $$$ 
% $$$ M = [];
% $$$ for i = 1:numP  
% $$$   A=[];
% $$$   for view = 1:numV
% $$$     PPM(:,:,view) = PPM(:,:,view)./norm(PPM(3,1:3,view));
% $$$     A =  [A
% $$$       PPM(1,:,view) - m(1,i,view)*PPM(3,:,view)
% $$$       PPM(2,:,view) - m(2,i,view)*PPM(3,:,view)];
% $$$   end   
% $$$   x = ns(A);
% $$$   M =[M x(1:3)./x(4)];
% $$$ end
% $$$ 
