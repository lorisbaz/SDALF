function [XX,YY] = draw_lines(L, X, Y, lt)
%DRAW_LINES Draw lines in cartesian form 
% draw_lines(lines)

% Author: Michael Bosse

if nargin < 4
   lt = '-';
end
if nargin < 3
   if nargin < 2 
      lt = '-';
   else
      lt = X;
   end
   V = axis;
   X = V(1:2); Y = V(3:4);

   hold on

end


i = find(abs(L(:,1)) < abs(L(:,2)));
j = find(abs(L(:,1)) >= abs(L(:,2)));

YY1 = -L(i,1)./L(i,2)*X + -L(i,3)./L(i,2)*ones(size(X));
XX1 = ones(length(i),1)*X;

XX2 = -L(j,2)./L(j,1)*Y + -L(j,3)./L(j,1)*ones(size(Y));
YY2 = ones(length(j),1)*Y;

XX = [XX1' XX2'];
YY = [YY1' YY2'];

[foo,k] = sort([i; j]);
plot(XX(:,k),YY(:,k),lt)

if nargin < 3
hold off
end
