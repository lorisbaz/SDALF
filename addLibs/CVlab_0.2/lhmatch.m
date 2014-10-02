function  match = lhmatch(d,sigma);
% LHMATCH Match two sets of points 
%% d e' la matrice delle distanze d(i,j)= distanza tra feature i e feature j
%% match(i) contiene l'indice della feature in corrispondenza (NaN altrimenti). 
%% match e' riferito allo spazio delle righe

%% per usarlo con distanza di mahalanobis (o comunque con distanze che gia'
%% tengono conto della varianza, porre sigma=1)

% d non deve nec. essere quadrata: se ci sono meno colonne che righe, ci sara'
% qualche feature nello spazio righe non assegnata (NaN nel vettore match);
% se ci sono meno colonne che righe ci sara' qualche feature nello spazio
% colonne non assegnata (non tutti gli indici colonna sono presenti in
% match).


% Algorithm: Scott & Longuett-HIggins 

% Author: Andrea Fusiello

match= NaN*ones(size(d,1),1);

G = exp(-d.^2/(2*sigma^2));
[U,S,V]  = svd(G);
P = U * eye(size(S)) * V';
% P is the same size as G

% Find dominant of each line and column
[V, J] = max(P');
[V, I] = max(P) ;

% 
% % Set a one-to-one correspondence between Ii and Jj
% % only if P(i,j) is MAX of both row i and column j

for i = 1:size(d,1)
  if I(J(i))==i
    match(i) = J(i);
  end 
end

