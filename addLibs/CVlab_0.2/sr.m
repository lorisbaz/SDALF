function [R,t] = sr(F,Al,Ar,ml,mr)
%SR calcola la fattorizzazione SR della matrice essenziale
%
%[R,t] = sr(F,Al,Ar,ml,mr) la matrice essenziale E e'
%calcolata trammite la matrice fondamentale F e le matrici dei
%parametri intrinseci Al, Ar. Poi viene fattorizzata in una matrice di
%rotazione R e un vettore di traslazione t.

%    Algorithm: Hartley, ECCV'92
%
%    Author: A. Fusiello 1999

%Controllo del formato dei parametri di input
[na,ma]=size(Al);
if na~=3 | ma~=3
    error('Formato errato della matrice dei parametri intrinseci (3x3)!!');
end

[rml,cml]=size(ml);
if (rml ~= 2)
    error('Le coordinate immagine devono essere cartesiane!!');
end

[rmr,cmr]=size(mr);
if (rmr ~= 2)
    error('Le coordinate immagine devono essere cartesiane!!');
end

[nF,mF]=size(F);
if nF~=3 | mF~=3
    error('Formato errato della matrice fondamentale (3x3)!!');
end


PPM = zeros(3,4,2);
m=zeros(2,size(ml,2),2);
m(:,:,1)=ml;
m(:,:,2)=mr;

% essential matrix
E = Ar'*F*Al;

%  SVD factorization such that E = U*D*V'
[U,D,V]= svd(E);

S1=[
    0 -1 0
    1 0 0
    0 0 0];

R1=[
    0 1 0
    -1 0 0
    0 0 1];

%  left perspective projection matrix
PPM(:,:,1) = Al *[
    1     0     0     0
    0     1     0     0
    0     0     1     0 ];

cheirality = false;
for j = 1:4

    %  skew symmetric matrix (representing translation)
    S= (-1)^j * U*S1*U';

    %  rotation matrix
    if j<=2
        R = det(U*V') * U*R1*V';
    end

    if j>2
        R = det(U*V') * U*R1'*V';
    end

    %  translation vector
    t=[S(3,2) S(1,3) S(2,1)]';
    t = t / norm(t);

    % 3D points distance from BOTH focal planes must be positive
    PPM(:,:,2) = Ar * [R   t];
    w=intersect(PPM,m,'midpoint');
    G0 = [R t];
    wt = rigid(G0,w);

    if min(w(3,:))>= 0 & min(wt(3,:))>= 0
        cheirality = true;
        break
    end

end

if ~cheirality
    warning('cheirality failed!')
end
