function [map_out_FG_intensity,map_out_FG_label,map_out_BG_intensity] =...
         Y_TestingImageMultiscaleMAG(img,F_arr,name,scale,sample,win_h,win_w,dimension)
% Y_TestingImageMultiscale testing LOGITBOOST Person detector on a single image.
%   PWatch performs person detection via covariance features on a
%   Riemaniann manifold.
%
%   [map_out_FG_intensity,map_out_FG_label,map_out_BG_intensity] =...
%         Y_TestingImageMultiscaleMAG(img,F_arr,name,scale,sample,win_h,win_w,dimension)
%
%   In:
%   img         := test image
%   F_arr       := classifier learnt (poroduced by Y_Training)
%   name        := path name for storing results
%   sample      := sample step on test image
%   win_h       := detection window size 
%   win_w       := detection window size 
%   dimension   := features number of covariance matrix
%   Out:
%   map_out_FG_intensity := scalar field (fg detection image)
%   map_out_FG_label     := scalar field  with per pixel class label of fg
%   map_out_BG_intensity := scalar field (bg detection image)

%   References:
%      [1]O. Tuzel, F. Porikli, P. Meer: Pedestrian Detection via
%      Classification on Riemannian Manifolds, Pattern Analysis and Machine
%      Intelligence, IEEE Transactions on : Accepted for future publication
%
%      [2] S. Kevin Zhou, J.H. Park, B. Georgescu, C. Simopoulos, J. Otsuki,
%      and D. Comaniciu: Image-based multiclass boosting and echocardiographic view
%      classification, IEEE Computer Society Conference on Computer Vision
%      and Pattern Recognition (CVPR 06)

%   Copyright 2009 Diego Tosato
%   $Revision: 1.3.6.1 $  $Date: 2009/08/28$

%% ---  General settings
d_sim            =   dimension*(dimension+1)/2+1;
img_or           =   img;
NumLevels        =   length(F_arr);
F                =   F_arr{1};
NumClass         =   length(F);
Ureg             =   [1,1,100,100];
l_size           =   zeros(1,NumClass); % number of weak learners per class

%% --- Setting detection wiondows
map_out_FG_intensity  =   zeros(size(img,1),size(img,2),length(scale),'single'); %result positive detection map
map_out_FG_label      =   zeros(size(img,1),size(img,2),length(scale),'single');%result positive detection map
map_out_BG_intensity  =   zeros(size(img,1),size(img,2),length(scale),'single');
maps_BG_intensity     =   zeros(size(img),'single'); %intensity result per level

% tic
%% --- Per scale analysis
for s = 1:length(scale)
    
    % --- scaling image (for multiscale analysis)
    %img_s    =   imresize(img_or,scale(s),'cubic');
    img            =   imresize(img_or,scale(s),'cubic');
    % --- padding
    [img,indd]    =   array_padd(img, [win_w/2, win_h/2], 0,'both','replicate'); %padding replicate


    maps_FG_intensity       =   zeros(size(img),'single'); %intensity result per level
    maps_FG_label           =   zeros(size(img),'single'); %label result per level
    feat                    =   Y_BuildF(img,[1,1,size(img,1),size(img,2)],dimension); % construct the map of features
    pxy                     =   zeros(size(feat)); % first order tensor of integral images
    Qxy                     =   zeros(size(feat,1),size(feat,2),size(feat,3),size(feat,3)); % second order tensor of integral images
    
%% --- Building covariance descriptor
    for ii = 1:dimension
        pxy(:,:,ii) =  cumsum(cumsum(feat(:,:,ii),2));
        for jj = ii:dimension
            Qxy(:,:,ii,jj)    =  cumsum(cumsum(feat(:,:,ii).*feat(:,:,jj)),2);
            Qxy(:,:,jj,ii)    =  Qxy(:,:,ii,jj);
        end
    end
    % analysis spatial limit
    range_h  =   size(feat,2) - win_h;
    range_w  =   size(feat,1) - win_w;

%% --- Main
    for j = 1:sample(s):range_h
        for i = 1:sample(s):range_w
            vote                =   zeros(1,NumClass);
%% --- covariance of C_R
            DataC       =  Y_CovCal_d(pxy,Qxy,[i,j,i+win_w-1,j+win_h-1],dimension,[i,j,i+win_w-1,j+win_h-1]);
            
            if any(diag(DataC) <= 0)
                D = diag(DataC);
                D(D <= 0) = 10^(-4);
                D = D - diag(DataC);
                DataC = DataC + diag(D);
            end
%% --- cascade analysis
            level = 1;
            while level <= NumLevels
                p                   =   ones(1,NumClass);
                goodClass           =   1:NumClass; % active node (= class) of TREE STRUCTURE
                p_merge             =   1; % merge probability (it is used in the tree structure)
                n_goodClass         =   length(goodClass);
                FkR(1:NumClass)     =   0; % final classification {class(sign)+confidence}
                F                   =   F_arr{level}; % load classifier at current level
                % --- number of wl per class
                for class = 1:NumClass
                    l_size(class)         =   size(F{class}.g,1);
                end
                % --- max number of wl
                l_size_max =  max(l_size);
%                 l_size_max = min(40,l_size_max);
                for ll=1:l_size_max
                    for class = 1:n_goodClass
%% --- coordinate calculation 
                        coords              =   F{goodClass(class)}.coords(ll,:);
                        r(1)                =   round(i + coords(1) * win_w/Ureg(3));
                        r(2)                =   round(j + coords(2) * win_h/Ureg(4));
                        r(3)                =   round(i + (coords(1) + coords(3)) * win_w/Ureg(3));
                        r(4)                =   round(j + (coords(2) + coords(4)) * win_h/Ureg(4));

                        % --- control
                        if r(1) < 1
                            r(1) = 1;
                        end
                        if r(2) < 1
                            r(2) = 1;
                        end
                        if r(3)-r(1)+1 >= win_w
                            r(3) = r(3) - (r(3)-r(1)+1 - win_w) -1;
                        end
                        if r(4)-r(2)+1 >= win_h
                            r(4) = r(4) - (r(4)-r(2)+1 - win_h) -1;
                        end
%% --- covariance of C_r
                        C                   = estimateCovarianceWindow(pxy,Qxy,r,dimension, DataC,[i,j,i+win_w-1,j+win_h-1]);
                        fk(goodClass(class))       =   [C',1]*F{goodClass(class)}.g(ll,:)';
                        %end
                    end

                    % --- control: dangerous!!!
                    %fk(fk < -sqrt(d_sim)) = -sqrt(d_sim);
                    %fk(fk > +sqrt(d_sim)) = +sqrt(d_sim);
                    % ---
                    
                    %% --- classifiers update
                    fk_sum = 1/(NumClass)*sum(fk(goodClass));
                    for class = 1:n_goodClass
                        fk(goodClass(class))        =   ((NumClass-1)/NumClass)*(fk(goodClass(class)) - fk_sum);
                        FkR(goodClass(class))       =   FkR(goodClass(class)) + fk(goodClass(class));
                    end
                    % --- posterior estimation
                    p_merge_l   =   0;
                    [p,p_merge_l] = X_ComputePosterior(p,n_goodClass,goodClass,FkR,p_merge,NumClass,p_merge_l,ll,l_size);
                    % ---  update p_merge
                    p_merge = p_merge .* (1-p_merge_l);

                    %--- delete class from tree structure
                    idx            = find(ll+1 > l_size(goodClass));
                    goodClass(idx) = [];
                    n_goodClass    = length(goodClass);
                end
                
                % --- multiclass classification
                [val_p,class_p] = max(p); % classification
                if class_p(1) ~= NumClass
                    vote(class_p) = vote(class_p) + 1;
                    max_vote    = max(vote);
                    parity      = sum(vote == max_vote);
                    if parity > 1
                        val         = val_p(1);
                        class       = class_p;
                    else
                        [val,class] = max(vote);
                        val         = p(class(1));
                    end
                else
                    val         = val_p(1);
                    class       = class_p;
                end

                II          = i+round(win_w/2); % detection on center of detection win
                JJ          = j+round(win_h/2);
                
                if class ~= NumClass
                    maps_FG_intensity(II-ceil(sample(s)/2):II+ceil(sample(s)/2),JJ-ceil(sample(s)/2):JJ+ceil(sample(s)/2)) =   val;
                    maps_FG_label(II-ceil(sample(s)/2):II+ceil(sample(s)/2),JJ-ceil(sample(s)/2):JJ+ceil(sample(s)/2))     =   class;
                else
                    maps_BG_intensity(II-ceil(sample(s)/2):II+ceil(sample(s)/2),JJ-ceil(sample(s)/2):JJ+ceil(sample(s)/2)) =  -val;
                    maps_FG_intensity(II-ceil(sample(s)/2):II+ceil(sample(s)/2),JJ-ceil(sample(s)/2):JJ+ceil(sample(s)/2)) =  0;
                    maps_FG_label(II-ceil(sample(s)/2):II+ceil(sample(s)/2),JJ-ceil(sample(s)/2):JJ+ceil(sample(s)/2))     =  0;
                    level = NumLevels; % not human -> exit
                end

                level = level +1;
            end
        end
    end
    
%% --- Partial result (only for k = K)
    %--- fg intensity
    maps_FG_tmp_intensity        =   maps_FG_intensity(indd(1):indd(2),indd(3):indd(4)); %delete padding
%     maps_FG_tmp_intensity        =   maps_FG_intensity; 
    maps_FG_tmp_intensity        =   imresize(maps_FG_tmp_intensity,[size(img_or,1) size(img_or,2)],'nearest'); % turn map on original scale
    map_out_FG_intensity(:,:,s)  =   maps_FG_tmp_intensity; % s <- scale; k <- level
    %--- fg label
    maps_FG_tmp_label            =   maps_FG_label(indd(1):indd(2),indd(3):indd(4)); %delete padding
%     maps_FG_tmp_label            =   maps_FG_label; 
    maps_FG_tmp_label            =   imresize(maps_FG_tmp_label,[size(img_or,1) size(img_or,2)],'nearest'); % turn map on original scale
    map_out_FG_label(:,:,s)      =   maps_FG_tmp_label; % s <- scale; k <- level
    
    %--- bg label
    maps_BG_tmp_intensity        =   maps_BG_intensity(indd(1):indd(2),indd(3):indd(4)); %delete padding
%     maps_BG_tmp_intensity        =   maps_BG_intensity;
    maps_BG_tmp_intensity        =   imresize(maps_BG_tmp_intensity,[size(img_or,1) size(img_or,2)],'nearest'); % turn map on original scale

    map_out_BG_intensity(:,:,s)  =   maps_BG_tmp_intensity; % s <- scale; k <- level
    
%     fprintf('Scale %1.2f done.. \n',scale(s));
    %----------------------------------------------------------
end


%% --- print result
% l_scale = length(scale);
% scrsz   = get(0,'ScreenSize');
% fig     = figure('Position',[1 scrsz(4)/2 scrsz(3)*2/3 scrsz(4)/3]);
% %subplot(1,l_scale+1,1);
% subplot(l_scale,3,1);
% imagesc(img_or); hold on, colormap(gray);axis equal, axis tight;title('Detection')%colorbar;
% %plot([1, win_h*(2-scale(i)), win_h*(2-scale(i))], [win_w*(2-scale(i)), win_w*(2-scale(i)), 1], 'b');
% map_out_mask = zeros(size(maps_FG_tmp_intensity(:,:,1)));
% for i=1:l_scale
%     Mpi = map_out_FG_intensity(:,:,i);
%     peak = max(max(Mpi));
%     if peak>0
%         [indx, indy] = find(Mpi > 0.0);
%         subplot(l_scale,3,i);
%         hold on
%         plot(indy, indx, 'rs',  'MarkerFaceColor','r','MarkerSize',0.5);axis equal, axis tight;
%          L        = Mpi;
%          % --- trashholding
%          %l        = graythresh(L);
%          %L(L>=l)   = 1; 
%          %L(L<l)    = 0; 
%          L = bwlabel(L);
%          stats    = regionprops(L, 'Centroid');
%          %centroid = round(stats(L).Centroid);
%          cx = [];cy = [];
%          for c = 1:size(stats,1)
%              cc = round(stats(c).Centroid);
%              [r_indx,r_indy] = find(L == c);
%              l_indx          = length(indx);
%              found           = 0;
%              l               = 1;
%              while l <= l_indx && found  == 0
%                  if  ~isempty(find(indx(l) == r_indx)) &&  ~isempty(find(indy(l) == r_indy))
%                      cx              = [cx;cc(2)];cy = [cy;cc(1)];
%                      found           = 1;
%                  end
%                  l = l + 1;
%              end
%          end
%         % --- 
%         for cc = 1:length(cx)
%             hold on
%             rectangle('Position',[cy(cc)-ceil(win_h*1/2),cx(cc)-ceil(win_w*1/2),win_h,win_w],'Curvature',[0,0],...
%                 'EdgeColor',[1 scale(i) scale(i)]);
%             if cx(cc)-ceil(win_w*1/2)< 1 &&  cy(cc)-ceil(win_h*1/2)< 1 
%                  map_out_mask(1:cx(cc)+ceil(win_w*1/2),1:cy(cc)+ceil(win_h*1/2)) = 1;
%             elseif cx(cc)-ceil(win_w*1/2)< 1
%                  map_out_mask(1:cx(cc)+ceil(win_w*1/2),cy(cc)-ceil(win_h*1/2):cy(cc)+ceil(win_h*1/2)) = 1;
%             elseif  cy(cc)-ceil(win_h*1/2)< 1   
%                  map_out_mask(cx(cc)-ceil(win_w*1/2):cx(cc)+ceil(win_w*1/2),1:cy(cc)+ceil(win_h*1/2)) = 1;
%             elseif cx(cc)+ ceil(win_w*1/2)>size(img_or,1) &&  cy(cc)+ceil(win_h*1/2)>size(img_or,2) 
%                  map_out_mask(cx(cc)-ceil(win_w*1/2):size(img_or,1),cy(cc)-ceil(win_h*1/2):size(img_or,2)) = 1; 
%             elseif cx(cc)+ ceil(win_w*1/2)>size(img_or,1)
%                  map_out_mask(cx(cc)-ceil(win_w*1/2):size(img_or,1),cy(cc)-ceil(win_h*1/2):cy(cc)+ceil(win_h*1/2)) = 1;
%             elseif cy(cc)+ceil(win_h*1/2)>size(img_or,2)
%                  map_out_mask(cx(cc)-ceil(win_w*1/2):cx(cc)+ceil(win_w*1/2),cy(cc)-ceil(win_h*1/2):size(img_or,2)) = 1;
%             else 
%                 map_out_mask(cx(cc)-ceil(win_w*1/2):cx(cc)+ceil(win_w*1/2),cy(cc)-ceil(win_h*1/2):cy(cc)+ceil(win_h*1/2)) = 1;
%             end
%         end
%         %name_1 = [name(1:end-5) '_c' '.mat']; %ATTENZIONE SI SALVANO I CENTROIDI QUI SOLO PER 1 SCALA DI ANALISI
%         %save(name_1,'cx','cy');      
%     end
%     %set(image(map_out_pos(:,:,s,level-1)), 'AlphaData', image(map_out_pos(:,:,s,level-1)));
%     subplot(l_scale,3,2);
%     imagesc(map_out_FG_intensity(:,:,i));axis equal, axis tight;title('Intensity map');colorbar;
%     subplot(l_scale,3,3);
%     imagesc(map_out_FG_label(:,:,i));axis equal, axis tight;title('Labels map');colorbar;
%     colormap jet;
%     %subplot(1,3,2+i);
%     %mesh(map_out_pos(:,:,i));
% end
% colormap gray
% Frame = getframe(fig);
% imwrite(Frame.cdata, name, 'jpeg');
% %name_1 = [name(1:end-5) '_fg' '.mat'];
% %save(name_1,'map_out_mask');
% close all;

% fprintf('\n');
% toc

function  CC = estimateCovarianceWindow(pxy,Qxy,r,d, DataC,R) 
% --- vector space settings
C = ones(d);
T = triu(C);
idx_triu = find(T ~= 0);
% --- covariance of Cr
%--- normalized covariance descriptor
C     =   diag(diag(DataC).^(-1/2))*Y_CovCal_d(pxy,Qxy,r,d,R)*diag(diag(DataC).^(-1/2))';
CC    =   C(idx_triu);

%--- singularities control
if any(isnan(CC))
    CC(isnan(CC)) = 0;
    disp('C <- NaN');
end

if any(isinf(CC))
    CC(isinf(CC)) = 0;
    disp('C <- Inf');
end


function  [p,p_merge_l]  = X_ComputePosterior(p,n_goodClass,goodClass,FkR,p_merge,J,p_merge_l,ll,l_size)
F_num                    =   exp(FkR(:,goodClass));
F_num (F_num==inf)       =   realmax;
F_den                    =   sum(exp(FkR(:,goodClass)),2);
F_den(F_den == inf)      =   realmax;
F_den                    =   repmat (F_den, 1, n_goodClass);
p_merge                  =   repmat (p_merge, 1, n_goodClass);

if n_goodClass == J+1
	p                     =   F_num./F_den; % probability of x being in class jth
else
	P_FG_tmp              =   F_num./F_den; % probability of x being in class jth
    p                     =   P_FG_tmp .* p_merge; % probability of x being in class jth
end

% --- update p_merge
p_merge_l   = p_merge_l + sum(p(l_size == ll));

p(p < 0) = 0;

