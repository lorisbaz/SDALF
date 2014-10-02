function G = exterior(A,model3d,data2d,method);
%EXTERIOR Solve exterior orientation
%    G = exterior(A,model3d,data2d) returns camera pose G given a list of
%    3D points (model3d), the corresponding 2D image points (data2d) and the
%    intrinsic parameters (A).
%
%    See also: RESECT

% Author: Andrea Fusiello



%Controllo del formato dei parametri di input
[nA,mA]=size(A);
if nA~=3 | mA~=3
    error('Formato errato della matrice dei parametri intrinseci (3x3)!!')
end

[rml,cml2]=size(data2d);
 
if (rml ~= 2)
     error('Le coordinate 3D  devono essere cartesiane!!');
 end


if (cml2<6)
     error('Ci devono essere almeno 6 punti corrispondenti!!');
    end
   
[rml,cml3]=size(model3d);

 if (rml ~= 3)
     error('Le coordinate immagine devono essere cartesiane!!');
 end
 
if (cml3<6)
     error('Ci devono essere almeno 6 punti corrispondenti!!');
    end
    
    
if (cml2 ~= cml3)
  error('i punti 2d e 3d devono essere in corrispondenza')
end

    

if nargin == 3
    method = 'fiore';
end

if strcmp(method,'fiore')
   G = exterior_fiore(A,model3d,data2d);
elseif strcmp(method,'dlt')
  G = exterior_dlt(A,model3d,data2d);
else
    error('metodo inesistente');
end

