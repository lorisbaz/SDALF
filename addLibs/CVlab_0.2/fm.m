function F = fm(m1,m2)
%FM Computes fundametal matrix from point correspondence
%
%F = fm(m1,m2) dati due insiemi di punti corrispondenti,
%m1 e m2, viene calcolata la matrice fondamentale F.
%I punti vengono standardizzati e viene utilizzato
%l'algoritmo lineare di Longuett-Higgins and Hartley.
%
%N.B.:
%
%Servono almeno 8 punti per il successo dell'algoritmo.

%    Algorithm: Hartley PAMI,1997
%
%    Author: A. Fusiello 1999


[rm1,cm1]=size(m1);
if (rm1 ~= 2)
    error('Le coordinate immagine devono essere cartesiane!!');
end

[rm2,cm2]=size(m2);
if (rm2 ~= 2)
    error('Le coordinate immagine devono essere cartesiane!!');
end

n_points = size(m1,2);
if n_points < 8
  error('Sono necessari almeno 8 punti!!!!!');
end

%pre-conditioning 
 
[SM1,m1] = precond2(m1);
[SM2,m2] = precond2(m2);

u1=m1(1,:)';
v1=m1(2,:)';
u2=m2(1,:)';
v2=m2(2,:)';

% Preparing the equation matrix; must have at least 9 rows
A = zeros(n_points,9);
A(:,1) = u2 .* u1;
A(:,2) = u1 .* v2;
A(:,3) = u1;
A(:,4) = v1 .* u2;
A(:,5) = v1 .* v2;
A(:,6) = v1;
A(:,7) = u2;
A(:,8) = v2;
A(:,9) = 1;

% solution vector corresponding to the 
% least singular value
f = ns(A);

% recover the matrix
F = reshape(f,3,3);

% apply the inverse scaling
F = SM2' * F * SM1;               

% Enforce singularity constraint on F
[U D V] = svd(F);
D(3,3) = 0;
F = U * D * V';

% normalize
F = F./norm(F);












