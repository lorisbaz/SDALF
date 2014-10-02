function [f_u,f_v,u,v,R,t]=parcam(P)
%PARCAM Extract parameters from camera matrix
%
%PARCAM Funzione che permette di estrarre i parametri intrinseci
%       ed estrinseci dalla matrice di proiezione prospettica.
%
%[R,t] = parcam(P) ritorna la matrice di rotazione (R) e il vettore di
%traslazione (t). (i parametri estrinseci).
%
%[f_u,f_v,u,v] = parcam(P) ritorna la distanza foncale orrizontale 
%espressa in pixel (f_u), la distanza focale verticale espressa in
%pixel (f_v), le coordinate del punto principale (u,v). (i parametri
%intrinseci).
%
%[f_u,f_v,u,v,R,t] = parcam(P) ritorna i parametri intrinsici ed estrinseci.
%
%N.B.
%
%1) Si suppone che  l'angolo tra gli assi u e v sia di 90° gradi.
%2) Si suppone che il piano immagine sia dietro il piano focale.

% Author: Andrea Fusiello

    %Controllo del formato dei parametri di input
    if (ismpp(P) == 0)
        warning('Matrice di proiezione prospettica non adeguata!!');
    end
    
    u=P(1,1:3)*P(3,1:3)'; 
    v=P(2,1:3)*P(3,1:3)'; 
    
    f_u=norm(cross(P(1,1:3),P(3,1:3)));
    f_v=norm(cross(P(2,1:3),P(3,1:3)));

    r3=P(3,1:3);
    r1=(P(1,1:3)-u*r3)*(1/f_u);
    r2=(P(2,1:3)-v*r3)*(1/f_v);
    
    t3=P(3,4);
    t1=(P(1,4)-u*t3)/f_u;
    t2=(P(2,4)-v*t3)/f_v;
    
    R=[r1;r2;r3];
    t=[t1;t2;t3];
    
    if nargout==2
        f_u=R;
        f_v=t;
    end
