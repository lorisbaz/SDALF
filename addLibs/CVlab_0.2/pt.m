function mt = pt(H,m);   
%PT applica una trasformazione proiettiva 2D o 3D.
%
%mt = pt(H,m) applica una omografia H del piano o dello spazio
%proiettivo a un insime di punti m.

% Author: Andrea Fusiello

if sum(size(H)) == 6
  % omografia del piano
  mt = p2t(H,m);
  elseif sum(size(H)) == 8
% omografia dello spazio
  mt = p3t(H,m);
else
  error('Trasformazione non implementata');
end

