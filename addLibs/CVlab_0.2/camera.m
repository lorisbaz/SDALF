function P=camera(A,lookp,eyep,up)
%CAMERA Build a camera matrix
%    P=camera(A,lookp,eyep,up) returns  a camera matrix P positioned at eyep,
%    pointing to lookp.  The  upward direction is "up"  and  A is the
%    intrinsic  parameters matrix.  In other terms, eyep is the optical
%    centre,  lookp is a point on the optical axis.  Note that up is a vector
%    that determines the orientation of the camera, but it is not equal to 
%    the y axis of the camera reference frame.
%
%    See also: RESECT

% Author: Andrea Fusiello, 2004


%P = camera(A,lookp,eyep,up) ritorna la MPP dati la matrice dei
%parametri intrinseci (A), il centro ottico (eyeP), un punto
%sull'asse ottico (lookp) e la direzione verticale (up).

    %Controllo del formato dei parametri di input
    [n,m]=size(A);
    if n~=3 | m~=3
        error('Formato errato della matrice dei parametri intrinseci (3x3)!!')
    end

    [rl,cl]=size(lookp);
    if rl~=3 | cl~=1
        error('Il lookp deve essere un vettore colonna di 3 elementi!!')
    end

    [re,ce]=size(eyep);
    if re~=3 | ce~=1
        error('Il eyep deve essere un vettore colonna di 3 elementi!!')
    end

    [ru,cu]=size(up);
    if ru~=3 | cu~=1
        error('Il up deve essere un vettore colonna di 3 elementi!!')
    end
    
    R(3,:) = ((lookp - eyep)/norm(lookp - eyep))'; %asse Z che passa tra il centro ottico e il look point
    R(1,:) = cross(up',R(3,:));
    R(1,:) = R(1,:)/norm(R(1,:)); %asse X ortogonale a Z e all'asse verticale
    R(2,:) = cross(R(3,:),R(1,:));%asse Y ortogonale al piano XZ

    
    P = A*[R, -R*eyep];    
