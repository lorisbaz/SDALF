function [S,mo] = precond(m)
%PRECOND Normalize the coordinates of a point set (2D or 3D). 
%
%   [S,mo] = PRECOND(m) scale and translate the 2D point coordinates so that 
%   the barycenter is in origin and the average radius  is sqrt(2)
%
%   See also ....

% Algorithm ref.: Hartley
%
% Author: A. Fusiello, 2006


[dimP,numP]=size(m);

if (dimP == 2)
  [S,mo] = precond2(m);
elseif (dimP == 3)
  [S,mo] = precond3(m);
  if (sum(m(3,:)) == numP)
    warning('Are you preconditioning 2D points?')
  end
else
  error('Wrong dimension');
end

   
