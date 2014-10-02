function [R,t] = horn(Al,Ar,ml,mr,R0,b)
%HORN Solve relative orientation with Horn's algorithm.
%
%[R,t] = horn(Al,Ar,ml,mr,R0,b) dati le matrice dei parametri intrinseci
%Al e Ar, un insieme di punti corrispondenti ml e mr, una matrice di
%rotazione iniziale RO (approssimazione della soluzione) e la 
%baseline b; calcola la trasformazione R,t.

% Author: Andrea Fusiello

%Controllo del formato dei parametri di input
[na,ma]=size(Al);
if na~=3 | ma~=3
    error('Formato errato della matrice dei parametri intrinseci (3x3)!!');
end
[na,ma]=size(Ar);
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

[nr0,mr0]=size(R0);
if nr0~=3 | mr0~=3
    error('Formato errato della matrice di rotazione (3x3)!!');
end

nb=length(b);
if nb~=3
   error('La baseline deve avere 3 elementi!!')
end


u1=ml(1,:)';
v1=ml(2,:)';
u2=mr(1,:)';
v2=mr(2,:)';

% constants
prevres = 1E40;
res = 1E39;
vanishing = 1E-16;

% if known, use here the variance on the determination of
% the vertical disparity of the points
sig0=1.0; % arbitrary
sig1=0.5; % vertical disparity variance for image 1
sig2=0.5; % vertical disparity variance for image 2

b = b/norm(b);
R = R0;

% change to normalized coordinates
mln = pt(inv(Al),[u1';v1']);
mrn = pt(inv(Ar),[u2';v2']);
u1=mln(1,:)';
v1=mln(2,:)';
u2=mrn(1,:)';
v2=mrn(2,:)';

% normalise ray vectors
r1 = zeros(3,size(u1,1));
for i=1:size(u1,1)
    r1(:,i) =  [u1(i) v1(i) 1]' / norm([u1(i) v1(i) 1]);
end

r2 = zeros(3,size(u2,1));
for i=1:size(u2,1)
    r2(:,i) =  [u2(i) v2(i) 1]' / norm([u2(i) v2(i) 1]);
end

% start iteration
iters=0;

while ( (prevres  - res  >  vanishing) )
    iters = iters+1;
    
    %  rotate left ray vectors
    r1p = (R * r1);
    
    % compute cross products
    c = zeros(3,size(u2,1));
    d = zeros(3,size(u2,1));
    t = zeros(1,size(u2,1));
    w =  diag(ones(1,size(u2,1)));
    
    for i=1:size(u2,1)
        c(:,i) = cross(r1p(:,i),r2(:,i));
        d(:,i) = cross(r1p(:,i), cross(r2(:,i), b ));
        t(1,i) = b' * c(:,i);  
        
        w(i,i) = (norm(c(:,i))^2 * sig0^2)/((cross(b,r2(:,i))'*c(:,i))^2 *...
            norm(r1p(:,i))^2*sig1^2 + (cross(b,r1p(:,i))'*c(:,i))^2 * norm(r2(:,i))^2*sig2^2);
        
    end
    
    
    % residual
    prevres=res;
    res = b'* ( c*c') *b;
    
    
    % compute solution
    L = [ (c*w)*c'   , (c*w)*d', b
        ((c*w)*d')', (d*w)*d', [0 0 0]'
        b'    , [0 0 0] , 0 ] ;
    
    y = - [(c*w)*t'
      (d*w)*t'
      0];
    
    a = L\y;
    
    % update baseline
    b=b+a(1:3);
    b = b/norm(b);
    
    % update rotation
    q = [ sqrt(1-0.25* norm(a(4:6))^2), 0.5*a(4:6)'];
    R = iquat(q(1),q(2:4)') * R;
    % enforce orthogonality of R
    [U,D,V] = svd(R);
    R = U * diag([1,1,1]) * V';
    
end

% results 
t = b;







