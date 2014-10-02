%% Detected head analisys: eliminate the heads that not correspond to the skin
load('head_ViPER');

reg             =   makecform('srgb2lab');
regback         =   makecform('lab2srgb');

addpath('../../Common LIB/')
addpath('../Symmetries/');

dir_a = 'D:/datasets/Tagging/VIPeR/cam_a';
dir_b = 'D:/datasets/Tagging/VIPeR/cam_b';

dir_a_list = dir(strcat(dir_a,'\*.bmp'));
dir_b_list = dir(strcat(dir_b,'\*.bmp'));

n_img       =   size(dir_a_list,1);
rgb_im      =   imread(strcat(dir_a,'\',dir_a_list(1).name));
[H,W,trash] =   size(rgb_im);

DIMW = 24;
h_range = [0.0 0.1]; 
s_range = [0.3 0.6]; 
search_range_H  =   [H/2-H/6,H/2+H/6];
search_range_W  =   [W/2-W/4,W/2+W/4];

hh = waitbar(0,'Image features computation...');
for i=1:n_img
%      figure; 
     for c=1:2
        if c==1; dire = dir_a; lista = dir_a_list; else; dire = dir_b; lista = dir_b_list; end
        img      = imread(strcat(dire,'/',lista(i).name));

        %% feature extraction
        img_hsv     =   rgb2hsv(img);
        tmp         =   img_hsv(:,:,3);
%         tmp         =   histeq(tmp); % Color Equalization
        img_hsv     =   cat(3,img_hsv(:,:,1),img_hsv(:,:,2),tmp); % eq. HSV
        img         =   hsv2rgb(img_hsv);
        img_cielab  =   applycform(img, reg); % eq. CIELAB
        
%         torso_legs = fminbnd(@(x) sym_dissimilar(x,img_cielab),1,H);
%         split_busto = fminbnd(@(x) sym_similar(x,img_cielab(1:torso_legs,:,:)),1,W);
%         split_ped  = fminbnd(@(x) sym_similar(x,img_cielab(torso_legs+1:end,:,:)),1,W);
		
		torso_legs = fminbnd(@(x) sym_dissimilar(x,img_cielab),search_range_H(1),search_range_H(2));
        split_busto = fminbnd(@(x) sym_similar(x,img_cielab(1:torso_legs,:,:)),search_range_W(1),search_range_W(2));
        split_ped  = fminbnd(@(x) sym_similar(x,img_cielab(torso_legs+1:end,:,:)),search_range_W(1),search_range_W(2));
        
        
% Smaller search space...
%         torso_legs = fminbnd(@(x) sym_dissimilar(x,img_cielab),search_range_H(1),search_range_H(2));
%         split_ped  = fminbnd(@(x) sym_similar(x,img_cielab),search_range_W(1),search_range_W(2));
        
         TLanti(i,c) = torso_legs;
         BUsim(i,c) = split_busto;
         LEGsim(i,c) = split_ped; 
         
         UP = img_hsv(1:TLanti(i,c),:,:);
         UPdx = UP(:,BUsim(i,c)+1:end,:);
         UPsx = UP(:,1:BUsim(i,c),:);
         ldx = W-BUsim(i,c)+1;  lsx = BUsim(i,c); stW = max(ldx,lsx);
         UPdxW = repmat(round(W/2):-1:round(W/2)-ldx+1,TLanti(i,c),1);
         UPsxW = repmat(fliplr(round(W/2):-1:round(W/2)-lsx+1),TLanti(i,c),1);
         UPW = [UPsxW,UPdxW];UPW =UPW./(max(UPW(:)));
         
         torso_head = fminbnd(@(x) wsym_dissimilar(x,UP,UPW),1,size(UP,1));
         HDanti(i,c) =torso_head;
         
         
     end

     waitbar(i/n_img,hh);
end
close(hh)

%% plot
maxplot = 24;
ii = 1;
for i=1:n_img
     for c=1:2
        if c==1; dire = dir_a; lista = dir_a_list; else; dire = dir_b; lista = dir_b_list; end
%         PLOTK  = [];
        img      = imread(strcat(dire,'/',lista(i).name));
%         subplot(maxplot*2/12,12,ii),imagesc(img), axis image
        if ~any(isnan(det_final(i,c,:)))
%             rectangle('Position',det_final(i,c,:),'EdgeColor',[0 1 0]);
%             text(det_final(i,c,1),det_final(i,c,2),'m','color',[0 0 0],'BackgroundColor',[1 1 1],'Margin',1,'FontSize',8);
        end
        UP_HD = img(1:HDanti(i,c),:,:); %imresize(rgb2gray(img(1:HDanti(i,c),:,:)),2,'nearest');
        map_out_FG_intensity = hs_detectSkin(UP_HD, h_range, s_range);
        BBox = regionprops(double(map_out_FG_intensity),'Centroid');
        if ~isempty(BBox)
            det_skin = [BBox.Centroid-DIMW/2 DIMW DIMW];
%             rectangle('Position',det_skin,'EdgeColor',[1 0 0]);
%             text(det_skin(1),det_skin(2),'s','color',[0 0 0],'BackgroundColor',[1 1 1],'Margin',1,'FontSize',8);
        else
            det_skin = [NaN,NaN,NaN,NaN];
        end
        if ~any(isnan(det_final(i,c,:))) && ~any(isnan(det_skin)) 
            % [TODO] intersection
            intersection = sum(rectint(det_final(i,c,:),det_skin));
            
            if intersection/(DIMW^2) >= 60/100 % detectors agree!
                det_final(i,c,1:2) = [(det_final(i,c,1)+det_skin(1))/2,(det_final(i,c,2)+det_skin(2))/2];
%                 rectangle('Position',det_final(i,c,:),'EdgeColor',[0 0 1]);
%                 text(det_final(i,c,1),det_final(i,c,2),'F','color',[0 0 0],'BackgroundColor',[1 1 1],'Margin',1,'FontSize',8);
            else
                det_final(i,c,:) = [NaN,NaN,NaN,NaN];
            end
        end
        ii = ii+1;
     end
%      if ii>=maxplot*2
%          ii = 1;
%          pause;
%      end
end