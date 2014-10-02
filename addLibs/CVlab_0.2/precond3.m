function [S,wo] = precond3(w)
%_PRECOND3 Normalize the coordinates of a 3D point set. 
%   
% Algorithm ref.: Hartley
%
% Author: A. Fusiello (adapted from original code by C. Plakas)


[rm2,cm2]=size(w);
if (rm2 ~= 3)
    error('Le coordinate  devono essere cartesiane!!');
end

x1=w(1,:)';
y1=w(2,:)';
z1=w(3,:)';


n_points = size(x1,1);


avgx1 = (sum(x1)/n_points);
avgy1 = (sum(y1)/n_points);
avgz1 = (sum(z1)/n_points);

x1 = x1 - avgx1;
y1 = y1 - avgy1;
z1 = z1 - avgz1;

tscale1 = sum(sqrt(x1.^2 + y1.^2 + z1.^2))/n_points/sqrt(3);

x1 = x1./tscale1;
y1 = y1./tscale1;
z1 = z1./tscale1;

S = [1/tscale1    0         0     -avgx1/tscale1; 
     0       1/tscale1      0     -avgy1/tscale1; 
     0            0    1/tscale1  -avgz1/tscale1;
     0            0        0              1];

wo = [x1';y1';z1'];

% end of scaling


% scale the point coordinates so that centroid is in origin 
% and average radius  is sqr(2)
% $$$ xc = (sum(M(1,:))/n_points);
% $$$ yc = (sum(M(2,:))/n_points);
% $$$ uc = (sum(m(1,:))/n_points);
% $$$ vc = (sum(m(2,:))/n_points);
% $$$ M(1,:) = M(1,:) - xc;
% $$$ M(2,:) = M(2,:) - yc;
% $$$ m(1,:) = m(1,:) - uc;
% $$$ m(2,:) = m(2,:) - vc;
% $$$ rM = sum(sqrt(M(1,:).^2 + M(2,:).^2))/n_points/sqrt(2);
% $$$ rm = sum(sqrt(m(1,:).^2 + m(2,:).^2))/n_points/sqrt(2);
% $$$ M(1,:) = M(1,:)/rM;
% $$$ M(2,:) = M(2,:)/rM;
% $$$ m(1,:) = m(1,:)/rm;
% $$$ m(2,:) = m(2,:)/rm;
% $$$ TM = [1 0 -xc; 
% $$$       0 1 -yc; 
% $$$       0 0  rM];
% $$$ Tm = [1 0 -uc; 
% $$$       0 1 -vc; 
% $$$       0 0  rm];
% $$$   
% $$$ % end of scaling
