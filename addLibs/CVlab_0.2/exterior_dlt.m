function G = exterior_dlt(A,model3d,data2d);
%_EXTERIOR_DLT solve exterior orientation with DLT
%    G = exterior(A,model3d,data2d) returns camera pose G given a list of
%    3D points (model3d), the corresponding 2D image points (data2d) and the
%    intrinsic parameters (A). Being an iterative algorithms, it requires a
%    starting guess (G0).
%
%    See also: RESECT

% Algorithm ref.: DLT

% Author: A. Fusiello, 2006


% change to normalized coordinates
% (or image coordinates) 
m = pt(inv(A),data2d);

G=resect(m,model3d);

% normalizzo, anche se lo faceva gia' calib
q3=G(3,1:3);
G=G./norm(q3);
    
% forzo la matrice di rotazione 
G(1:3,1:3) = plr(G(1:3,1:3)) ;      
    
    
    
    





