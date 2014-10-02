% MS_CLUSTERING does Mean Shift based clustering.
%
%   [Label,C_centers,y]= MS_Clustering(Data,h) clusters the data point in Data(num_data,num_dim)
%   using a Mask_type kernel mask, equal for each dimension. The width of the mask 
%   (possible different) for each dimension is located in h.
%   - The result is given in Labels(num_data,1), in which each data point is labeled 
%   with the number c of the correspondent cluster. 
%   - C_Centers represents the
%   center for each c cluster in the num_dim-th dimensional space
%   In the current version Mask_type has to be 'Uniform'.  
%   Applico un size filter, regolato da min_size che da la numerosita'
%   minima di un cluster
function [Label,C_centers,nearest_data]= MS_Clustering2(Data,h,Mask_type,subspace,min_size)
if ~strcmp(Mask_type,'Uniform')
    fprintf('uncorrect kernel shape (in Mask_type)!!!');
    Label                  =   [];
    return
end
[num_data, total_dims]     =   size(Data); 
num_subspace          =   length(subspace);

% --------------------
% MS_CLUSTERING
% --------------------
%hw                      =   waitbar(0,'performing MS Clustering ...');
Right_ones              =   ones(num_data,1);
for n=1:num_data
    %waitbar(n/num_data,hw)
    yNow                    =   Data(n,:);
    delta                   =   1000000; 
    t                       =   0;
    iteration =0;
    while delta>0.000000001
        if iteration ==100000; break;end
        t                       =   t+1;
        beg                     =   0;
        for sub=1:num_subspace % in num_subspace ho il numero di sottospazi in cui mi muovo
            num_dim = subspace(sub); % il numero di dimensioni del sub-esimo sottospazio
            for dim=(beg+1):(beg+num_dim)
                Ris(:,dim)              =   (Right_ones.*yNow(dim))-Data(:,dim);
                Ris(:,dim)              =   Ris(:,dim)./h(dim);
                Ris(:,dim)              =   Ris(:,dim).^2;
            end
            Evalu                   =   Ris(:,(beg+1):(beg+num_dim));
            Norma                   =   sum(Evalu,2);
            Sum_KerG(:,sub)         =   Norma<=1;
            beg                     =   beg + num_dim;
        end
            Kernel_response         =   repmat(prod(double(Sum_KerG),2),1,total_dims);
            Weight                  =   sum(Data.*Kernel_response,1); 
            Weight_den              =   sum(Kernel_response,1); 
            yNext                   =   Weight./Weight_den;
            %yNext(1:3)              =   yNext(1:3)/norm(yNext(1:3));     
            delta                   =   sqrt(sum(abs(yNext-yNow).^2));
            yNow                    =   yNext;
            iteration               =   iteration + 1;
    end
    End_data_position(n,:)  =   yNow;
end
%close(hw)

% --------------------
% DATA LABELING
% --------------------

center_radius   =   0.8*h;
Unknow          =   1:num_data; % I punti ancora da assegnare
clust_index     =   1;
Label           =   zeros(num_data,1);

while ~isempty(Unknow)
    num_Unknow  =   length(Unknow); % All'inizio sono tutti i punti;
    n           =   Unknow(1); % Il centro che vado a valutare e' il primo dei punti
                               % ancora senza etichetta;
    Temp_center =   End_data_position(n,:); % Ricavo le ccordinate di tale centro;
    To_process  =   End_data_position(Unknow,:); % Prendo le coordinate di quelli ancora da valutare;
    Ris         =   zeros(num_Unknow,total_dims); % Creo un vettore solo per i punti in analisi.
    beg         =   0; % indice con cui indico i vari sottospazi che sto analizzando;
    Good        =   zeros(num_Unknow,num_subspace);
    for sub=1:num_subspace % in num_subspace ho il numero di sottospazi in cui mi muovo
        num_dim = subspace(sub); % il numero di dimensioni del sub-esimo sottospazio
        for dim=(beg+1):(beg+num_dim)
            Ris(:,dim)= ((ones(num_Unknow,1).*Temp_center(dim))-To_process(:,dim)).^2;
        end
        Evalu                   =   Ris(:,(beg+1):(beg+num_dim));
        Norma                   =   sqrt(sum(Evalu,2));
        Good(:,sub)             =   abs(Norma)<=center_radius(dim);
        beg                     =   beg + num_dim;
    end
    Proximity_response  =   prod(double(Good),2);
    Neigh               =   find(Proximity_response==1);
    if length(Neigh)>min_size  %-----size filter
        Label(Unknow(Neigh))=   clust_index;
        clust_index         =   clust_index + 1;
    end
    Unknow(Neigh)       =   0;
    Unknow              =   Unknow(find(Unknow~=0));
end

% --------------------
% CENTERS COMPUTATION
% --------------------
    C_centers=zeros(clust_index-1,total_dims);
for c=1:clust_index-1
     Good_index = find(Label==c);
     C_centers(c,:) = sum(End_data_position(Good_index,:),1)./length(Good_index);
end

for c=1:clust_index-1
    beg=0;
    Good_index = find(Label==c);
    dist_par = [];
     for sub=1:num_subspace % in num_subspace ho il numero di sottospazi in cui mi muovo
            num_dim = subspace(sub); % il numero di dimensioni del sub-esimo sottospazio
            for dim=(beg+1):(beg+num_dim)
                dist_par(:,sub)=sqrt(sum((repmat(C_centers(c,dim),length(Good_index),1)-Data(Good_index,dim)).^2,2));
            end
            beg                     =   beg + num_dim;
     end
     dist_clus =prod(dist_par,2);
     [min_dist,nearest] = min(dist_clus);
     nearest_data(c)=Good_index(nearest);
end
    
      
        
    
    
    