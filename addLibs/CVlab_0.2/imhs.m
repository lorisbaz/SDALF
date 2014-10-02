function [X,Y] = imimhs(I,ws,MAXFRAC,NP)
%IMHS Extracts salient points with Harris&Stephens' algorithm 
%   
%
%[X,Y] = imhs(I,ws,MAXFRAC,NP)
%MAXFRAC e' la soglia rispetto al massimo
%NP e' il numero di punti che si vogliono
%se NP=0 vengono restituiti tutti quelli
%che superano MAXFRAC
%I e' l'immagine e ws è la dimensione della
%matrice del filtro gaussiano.


% Author: Andrea Fusiello
% 
calcola derivate direzionali
Sob = fspecial('sobel'); 
Iu = filter2(Sob, I, 'same');
Iv = filter2(Sob',I, 'same');

% calcola convoluzione con gaussiana
% e le componenti della matrice 2x2 C
G = fspecial('gaussian',[ws,ws]); 
Iuv = filter2(G, Iu.*Iv, 'same');
Ivv = filter2(G, Iv.^2, 'same');
Iuu = filter2(G, Iu.^2, 'same');

% traccia e determinante
TR = Iuu + Ivv;
DET = Iuu.*Ivv - Iuv.^2;

% noble
%%C = (DET)./(TR.^2);
% harris&stephens
C = DET - 0.04 *TR.^2;

% sopprimo bordi
C(1:ceil(ws/2),:) =0;
C(end-ceil(ws/2):end,:) =0;
C(:,1:ceil(ws/2)) =0;
C(:,end-ceil(ws/2):end) =0;

% sopprime massimi locali
A=colfilt(C,[5 5],'sliding',@max);
local_maxima = find(C<A);
C(local_maxima) = 0;

% trova i valori piu' alti 
m = max(C(:));
i = find(C>MAXFRAC*m);  

if NP > 0
    if NP>length(i) 
        warning('non ci sono abbastanza punti') 
        NP=length(i);
    end
    %Sort the output and pick the top ones
    [foo,ind] = sort(-C(i));
    k = i(ind(1:NP));  
else
    k=i;
end

%Convert linear index into coordinates of the peaks.
[X,Y] = ind2sub(size(C),k);  
