function piramide(m)

% dati 4 punti 2d o 3d li disegna assumendo che siano i vertici di un tetraedro

% Author: Andrea Fusiello

if size(m,1) == 3
    
    x=m(1,:);
    y=m(2,:);
    z=m(3,:);
    
    x=[x(1) x(2) x(3) x(1) x(4) x(2) x(3) x(4)];
    y=[y(1) y(2) y(3) y(1) y(4) y(2) y(3) y(4)];
    z=[z(1) z(2) z(3) z(1) z(4) z(2) z(3) z(4)];
    
    plot3(x,y,z);
    
elseif size(m,1) == 2
        
        x=m(1,:);
        y=m(2,:);
        x=[x(1) x(2) x(3) x(1) x(4) x(2) x(3) x(4)];
        y=[y(1) y(2) y(3) y(1) y(4) y(2) y(3) y(4)];
   
        plot(x,y);
         hold on
        plot([x(1) x(2) ],[y(1) y(2) ],'-y')
        plot([ x(2) x(3) ],[ y(2) y(3) ],'-m')
        plot([ x(3) x(1)],[ y(3) y(1)],'-c')
        hold off
        
else
    error('Non supportato');
end
