function r = mrmse(PPM,m,w)
%MRMSE Multiple view RMSE
%
%   r = MRMSE(P,m,w) computes the root mean square error   
%   between image points m and the projection of the corresponding
%   3D points w with camera matrices in PPM.
%    
%
%   See also RMSE

% Algorithm ref.: none
%
% Author: A. Fusiello, 2006

siz = size(m);

if (size(siz,2) ~= 3)
error('Formato errato per i dati')
end

numberOfViews = siz(3);

r=0;

for view = 1:numberOfViews
    mp = proj(PPM(:,:,view),w);
    r = r + rmse(mp,m(:,:,view)) ;
end

r = r / sqrt(numberOfViews);
