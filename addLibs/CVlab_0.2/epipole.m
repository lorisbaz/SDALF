function e = epipole(F)
%EPIPOLE Extract the left epipole from a fundamental matrix
%    e = epipole(F) return the epipole corresponding to F such that F*e=0
%
%    See also: FUND


% Author: Andrea Fusiello


    %Controllo del formato dei parametri di input
    [n,m]=size(F);
    if n~=3 | m~=3
        error('Formato errato della matrice fondamentale (3x3)!!')
    end

    [U,D,V]= svd(F);

    if abs(V(3,3)) < 10^(-8)
        warning('epipolo all''infinito!!!');
    end

    e = V(:,3);
    
    e = e./norm(e);
