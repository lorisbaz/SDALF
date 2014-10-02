function enlEpitome = enlargeEpitome( epitome, fac )
% enlEpitome = enlargeEpitome( epitome, fac )
% 
% epitome : mean or variance of the epitome
% fac [0 1]: percentuale dell'epitomo da aggiungere ad ogni parte 
%
% e.g. muEnlarged = enlargeEpitome( mu, .2 )
% Zonta il 20% 

[Re Ce foo] = size(epitome);
aRe = ceil( fac.*Re );
aCe = ceil( fac.*Ce );
Up = epitome(1:aRe,:,:);
Lo = epitome(end-aRe+1:end,:,:);
Le = epitome(:,1:aCe,:);
Ri = epitome(:,end-aCe+1:end,:);

epitome3 = cat(2,Ri,epitome,Le);
bandU = cat(2,zeros(size(Up,1),size(Ri,2),3),Up,zeros(size(Up,1),size(Le,2),3));
bandL = cat(2,zeros(size(Up,1),size(Ri,2),3),Lo,zeros(size(Up,1),size(Le,2),3));
enlEpitome = cat(1,bandL,epitome3,bandU);
