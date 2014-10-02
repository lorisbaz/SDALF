function [PPM, w] = prec(m)
%PREC compute a projective reconstruction (structure and motion)
%     da un insieme di punti immagine corrispondenti.
%

% Algorithm: Sturm & Triggs (multiview)

% Author: Andrea Fusiello

siz = size(m);

numP = siz(2);
numV = siz(3);


if numV==2
  % ricostruzione da due viste

  F = fm(m(:,:,1),m(:,:,2));

  e = epipole(F');

  PPM(:,:,1) = [eye(3), [0 0 0]'];
  PPM(:,:,2) = [star(e)*F,  e];

  % sanity check
  % $$$     [FF,aa,bb] = fund(PPM(:,:,1),PPM(:,:,2));
  % $$$     s = svd([F(:),FF(:)]);
  % $$$     fprintf('Errore algebrico: %0.5g \n',s(2) );


  % midpoint va bene solo per ric. metrica
  w=intersect(PPM,m,'le');


else
  %  ricostruzione da n viste
   
  [PPM, w] = prec_triggs(m);
  
  %         ZULIANI 
  % $$$       
  % $$$      visibility = ones(numV,numP);
  % $$$   
  % $$$      [P, S, iter, E, x_star] = projective_reconstruction(m, visibility, ...
  % $$$ 	 1e-16, 200, 1e-16, 1); 
  % $$$   
  % $$$      % converto in coord. cartesiane
  % $$$      w = S(1:3,:)./ [S(4,:); S(4,:); S(4,:)];
  % $$$    
  % $$$      % metto le PPM nel modo giusto
  % $$$      for j = 1:numV
  % $$$        PPM(:,:,j) =  P(3*j-2:3*j,:);  
  % $$$      end

end

return


% $$$ function Mn = col_norm(M)
% $$$ 
% $$$ % normalizzo le colonne   
% $$$ Mn = [];
% $$$ 
% $$$ for i = 1:size(M,2)
% $$$ 
% $$$   Mn = [Mn, M(:,i)./norm( M(:,i) )];
% $$$ 
% $$$ end
% $$$ 
% $$$ return	 
% $$$ 
% $$$ function Mn = row_norm(M)
% $$$ 
% $$$ % normalizzo le righe
% $$$ Mn = [];
% $$$ 
% $$$ for j = 1:size(M,1)
% $$$ 
% $$$   Mn = [Mn; M(j,:)./norm( M(j,:) )];
% $$$ 
% $$$ end
% $$$ 
% $$$ return



function [PPM, w] = prec_triggs(m)

siz = size(m);
numP = siz(2);
numV = siz(3);


% preconditioning (it is crucial!)
for j = 1:numV
  [SM(:,:,j),m(:,:,j)] = precond(m(:,:,j));
end


% build the measurement matrix
M1 = [];
for i = 1:numP
  a = [];
  for j = 1:numV
    a = [a; m(:,i,j);1];
  end
  M1 = [M1, a];
end

% $$$    DD = [];
% $$$    for j = 1:numV
% $$$         D = [];
% $$$            for i = 1:numP
% $$$                D =  [D
% $$$                    zeros(3,i-1) [m(:,i,j);1] zeros(3,numP-i)];
% $$$            end
% $$$        DD(:,:,j) = pinv(D);
% $$$   
% $$$            end

% auxiliary matrix containing measurements 	 
QQ = [];
iQQ = [];
for i = 1:numP
  % Q_i contiene le proiezioni dell'i-esimo punto nelle varie viste
  Q_i = [];
  for j= 1:numV
    Q_i =  [Q_i
      zeros(3,j-1) [m(:,i,j);1] zeros(3,numV-j)];
  end
  QQ(:,:,i) = Q_i;
  iQQ(:,:,i) = pinv(Q_i);
end


niter = 0;
res = 1e40;
M = M1; % inizializza misure con depth=1


while (niter < 100)
  
  % the measurement matrix must be normalized to avoid trivial solutions
  M ./ norm(M,'fro');
    
  % trova struttura e moto date le misure M
  [U,D,V] = svd(M);
  P = U(:,1:4)*D(1:4,1:4);
  S = V(:,1:4)';

  % stima le depth z e aggiorna le misure

  % $$$          M = [];
  % $$$ 	 % stima per righe
  % $$$          for j = 1:numV
  % $$$               z = DD(:,:,j)* vec(P(3*j-2:3*j,:) * S);
  % $$$               M = [M
  % $$$                  M1(3*j-2:3*j,:) * diag(z)];
  % $$$           end

  M = [];
  % stima per colonne 
  for i = 1:numP
    z = iQQ(:,:,i)* P * S(:,i);	     
    M = [M, QQ(:,:,i)*z;];
  end
     
  % aggiornamento del residuo e controllo convergenza
  niter = niter + 1;
  prevres = res;
  res = norm(M-P*S,'fro');
  if abs(res - prevres) < eps
    break
  end
  if res > prevres
    warning('increasing residual');
  end
    
end

% normalizza i punti 3D w 
w = S(1:3,:)./ [S(4,:); S(4,:); S(4,:)];

% mette le matrici nella forma richiesta
for j = 1:numV
  PPM(:,:,j) = P(3*j-2:3*j,:);
end

%% fprintf('Errore di ricostruzione nelle immagini: %0.5g pixel \n',mrmse(PPM,m,w));

% postconditioning
for j = 1:numV
  PPM(:,:,j) = inv(SM(:,:,j)) * PPM(:,:,j);
end

% normalize such that first camera is = [I|0]
G0 = [PPM(:,:,1); [ 0 0 0 1]];
% apply  to cameras ...
iG0 = inv(G0);
for j = 1:numV
  PPM(:,:,j) =  PPM(:,:,j) * iG0;
end 
% ... and to points as well
w = p3t(G0,w);

% fine 
return
