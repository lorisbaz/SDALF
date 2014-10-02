function imo=imstereo(ims,imd,drange,ws)
%IMSTEREO Stereo block-matching with sum of squared differences
%
%imo = imstereo(nomeIs,nomeId,drange,ws) calcola una mappa di disparita'
%imo, con eventuali NaN dove non e' definita. I parametri sono:
%nomeIs: nome file dell'immagine sinistra
%nomeId: nome file dell'immagine destra
%drange:  intervallo di ricerca disparita' [di,ds]
%ws: dimensioni finestra di correlazione [a,b].

% Author: Andrea Fusiello

if (length(drange)~=2)
    error('L''intervallo di disparita'' deve avere 2 valori!!')
end
if (length(ws)~=2)
    error('La finestra di correlazione deve avere 2 valori!!')
end

dmin=drange(1);
dmax=drange(2);

ave=fspecial('average',ws);
[righe,colonne]=size(ims);
c=ones(righe,colonne,dmax-dmin+1).*Inf;

for d=0:dmax-dmin
    pad=ones(righe,dmin+d).*NaN;
    imst=[pad,ims(:,1:colonne-(dmin+d))]; %immagine traslata
    diff=(imd-imst).^2;
    c(:,:,d+1)=filter2(ave,diff);
    %imshow(diff,[]);
    %pause;   
end; 
[val,imo]=min(c,[],3);


imo(find(isnan(val)))=NaN;  
imo=imo+dmin-1;           
