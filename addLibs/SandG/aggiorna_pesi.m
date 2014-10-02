function [wR , muR, p_stdR] = aggiorna_pesi(x,w,mu,p_std, c_matched, alfa,numComps,min_p_std,max_p_std)

% x = [1,2,3]';
% w=[10,11,12;20,21,22;30,31,32];
% mu=w;
% p_std =[1,0,0;2,1,0;3,0,0];
% c_matched = [1;2;1];
% numComps = 3;

X=size(w,1); % Numero di osservazioni; 
Y=numComps;
x=x*ones(1,numComps);% ripeto le osservazioni per il numero di componenti;
first=X*Y; % 

w=reshape(w,first,1);
mu=reshape(mu,first,1);
p_std=reshape(p_std,first,1);
x=reshape(x,first,1); % Allineo tutti i valori per un unica colonna;

s_matched=Couple2single((1:X)',c_matched,X,Y); 
% Trovo le posizioni da modificare all'interno di questo vettore;
w=w.*(1-alfa);
w(s_matched)=w(s_matched)+ alfa;
diff=(x(s_matched)-mu(s_matched)).^2;
ro=alfa.*(1./sqrt(2*pi*p_std(s_matched))).*exp(-0.5.*diff./(p_std(s_matched).^2));
% if ro ==0
%     fprintf('niente!!!')
% end
mu(s_matched)=mu(s_matched).*(1-ro)+ro.*x(s_matched);
p_std(s_matched)=p_std(s_matched).*(1-ro)+ro.*(diff);% controllo sulla varianza minima !!!
p_std(s_matched) = (p_std(s_matched)<min_p_std)*min_p_std + (p_std(s_matched)>=min_p_std).*p_std(s_matched);
p_std(s_matched) = (p_std(s_matched)>max_p_std)*max_p_std + (p_std(s_matched)<=max_p_std).*p_std(s_matched);
wR=reshape(w,X,Y);
p_stdR=reshape(p_std,X,Y);

muR=reshape(mu,X,Y);
Sum=sum(wR');
Sum2=(Sum'*ones(1,numComps));
Sum2=reshape(Sum2,first,1);%[Sum';Sum';Sum'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AG1: schema di aggiornamento 1
w=w./Sum2; % Aggiornamento simultaneo. Dopo la prima iterazione, 
% le gaussiane centrate sui valori di pixel
% osservati hanno peso 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AG2: schema di aggiornamento 2
%  good = Sum2>1;
%  w(good)=w(good)./Sum2(good); 
% Cosi' faccio in modo che all'inizio non avvenga un'inizializzazione
% con solo un'iterazione (i pesi della componente gaussiana che matchano la
% prima iterazione diventano 1 la seconda, mentre io voglio che
% l'importanza cresca gradualmente;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wR=reshape(w,X,Y);
