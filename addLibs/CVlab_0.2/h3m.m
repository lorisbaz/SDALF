function  H=h3m(w,wp)
%_H3M Compute 3D homography btw two point sets
%
%H = h3m(w,wp) dati i punti mondo w, e i punti mondo 
%trasformati wp, calcola la matrice di traformazione H
%  wp = H*w

% Author: Andrea Fusiello

[rw,cw]=size(w);
if (rw ~= 3)
    error('Le coordinate mondo devono essere cartesiane!!');
end

[rwp,cwp]=size(wp);
if (rwp ~= 3)
    error('Le coordinate mondo devono essere cartesiane!!');
end


dime = size(w,2);
x=w(1,:)';
y=w(2,:)';
z=w(3,:)';

xp=wp(1,:)';
yp=wp(2,:)';
zp=wp(3,:)';

%costruisco la matrice del sistema

A =[];
b =[];
for i = 1:dime
  col = [ 
      x(i)        ,0           ,0
      y(i)        ,0           ,0
      z(i)        ,0           ,0
      1	          ,0           ,0
      0	          ,x(i)        ,0
      0	          ,y(i)        ,0
      0           ,z(i)        ,0
      0	          ,1           ,0
      0           ,0           ,x(i)
      0           ,0           ,y(i)
      0           ,0           ,z(i)
      0           ,0           ,1
     -x(i)*xp(i)  ,-x(i)*yp(i) ,-x(i)*zp(i)
     -y(i)*xp(i)  ,-y(i)*yp(i) ,-y(i)*zp(i)
     -z(i)*xp(i)  ,-z(i)*yp(i) ,-z(i)*zp(i)];
 
 noto = [xp(i), yp(i), zp(i)];
 
  A = [A col];
  b = [b noto];
  
end
A=A';
b=b';

q = pinv(A) * b;


H(1,1)=q(1) ;
H(1,2)=q(2) ;
H(1,3)=q(3) ;
H(1,4)=q(4) ;

H(2,1)=q(5) ;
H(2,2)=q(6) ;
H(2,3)=q(7) ;
H(2,4)=q(8) ;

H(3,1)=q(9) ;
H(3,2)=q(10) ;
H(3,3)=q(11) ;
H(3,4)=q(12) ;

H(4,1)=q(13) ;
H(4,2)=q(14) ;
H(4,3)=q(15) ;
H(4,4)=1 ;
