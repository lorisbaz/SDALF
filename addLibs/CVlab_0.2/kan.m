function R = kan(r)
%KAN calcola la matrice di rotazione di un certo angolo
%     intorno ad un asse dello spazio.
%
%R = kan(r) ritorna la matrice di rotazione di un angolo
%phi rispetto all'asse caratterizzato dal versore r.

% Author: Andrea Fusiello

nr=length(r);
if (nr ~= 3)
    error('Il versore r deve essere un vettore di 3 elementi!!');
end

phi = norm(r);
r = r/phi;

c = cos(phi);
s = sin(phi);
v = 1-cos(phi);

r=r/norm(r);

R = [ 	r(1)*r(1)*v+c,		r(1)*r(2)*v-r(3)*s,	r(1)*r(3)*v+r(2)*s
	r(1)*r(2)*v+r(3)*s,	r(2)*r(2)*v+c,		r(2)*r(3)*v-r(1)*s
	r(1)*r(3)*v-r(2)*s,	r(2)*r(3)*v+r(1)*s,	r(3)*r(3)*v+c  ];
