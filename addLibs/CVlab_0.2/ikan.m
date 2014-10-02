function [r] = ikan(R)
%IKAN Computes the axis and angle of a rotation matrix
%
%IKAN  individua l'asse e l'angolo di rotazione connessi ad
%      una generica matrice di rotazione.
%
%[r,phi] = ikan(R) ritorna il versore r che caratterizza l'asse
%e l'angolo phi connessi ad una generica matrice di rotazione.
%
%N.B.:
%Per angoli di rotazioni molto piccoli, prossimi allo zero,
%si ha una singolarita' questo problema viene risolto con i
%quaternioni.

% Author: Andrea Fusiello

%Controllo che sia una matrice di rotazione
if (~isrot(R))
    error('L''argomento non e' una matrice di rotazione!!')
end

phi = acos((trace(R)-1)/2);
r0 = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]';
r = phi* r0/norm(r0);
