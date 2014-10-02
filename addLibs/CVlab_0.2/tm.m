function T = tm(m1,m2,m3)
%TM calcola la matrice trifocale da corrispondenze
%
% Author: Andrea Fusiello

numP = size(m1,2);

% con il pre-conditioning  migliora il condizionamento di A, ma peggiore
% l'errore algeberico della ricostruzione ... che sia sbagliato il post?

% pre-conditioning 
% $$$ [H1,m1] = precond(m1);
% $$$ [H2,m2] = precond(m2);
% $$$ [H3,m3] = precond(m3);

A = [];
for i = 1:numP
  A =  [A
    kron(kron([m1(:,i);1]', star([m3(:,i);1])),star([m2(:,i);1]))];
end

t = ns(A);
T = vtrans(t,9);

% post-conditioning
% T = kron(inv(H3), inv(H2)) * T * H1;

