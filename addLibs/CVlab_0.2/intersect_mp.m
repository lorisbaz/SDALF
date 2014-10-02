function wr = intersect_mp(PPM,imagePoints)
%_INTERSECT_MP Intersect with midpoint algorithm
%
%    Given a set of projection matrices (PPM) and a set of corresponding
%    image points (the position gives the correspondence) the 3D position of
%    the 3D points that projects onto the image points is computed. Since
%    rays may not intersect, the mid-point of the common perpendicular to
%    the two rays is chosen (the mid-point method)
%
%    Do not use with projective reconstruction.


%    Algorithm: Beardsley, Ziserman, Murray IJCV 23(3),1997
%
%    Author: A. Fusiello 1999

%%fprintf('MP\n');

siz = size(imagePoints);
numberOfPoints = siz(2);
numberOfViews =  siz(3);
struc = [];
I = eye(3,3);

for i = 1:numberOfPoints
    
    A =  zeros(3,3);
    b =  zeros(3,1);
    
    for view = 1:numberOfViews
        
        Q =  inv(PPM(:,1:3,view));
        c =  -Q*PPM(:,4,view);
        d =   Q*[imagePoints(:,i,view); 1];
        d  = d/norm(d);
        
        A = A + I - d*d';
        b = b + c - (c'*d)*d;
        
    end
    
    w = inv(A) * b;
    struc = [struc  w];
   
end


wr=struc;






