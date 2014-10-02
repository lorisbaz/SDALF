function h =  plotref(G,s,color)
%PLOTREF plot a reference frame with the specified color and scale

% Author: Andrea Fusiello

if nargin == 1
    color = 'r';
    s = 1;
end

if nargin == 2
    color = 'r';
end


R = G(1:3,1:3);
t = G(1:3,4);
eyep = (-R'*t)';

a = (eyep);
plot3(a(1),a(2),a(3),['o', color]);
b = [s*R(1,:)+eyep;
     s*R(2,:)+eyep;
     s*R(3,:)+eyep];
plot3([a(1), b(3,1)],[a(2), b(3,2)],[a(3), b(3,3)],['-', color]);
plot3([a(1), b(2,1)],[a(2), b(2,2)],[a(3), b(2,3)],['-', color]);
h = plot3([a(1), b(1,1)],[a(2), b(1,2)],[a(3), b(1,3)],['-', color]);


