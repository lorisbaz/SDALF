%% UPLOADING AND FEATURE EXTRACTION
% MSK = ones(128,64); % COMMENT THIS!!!! (for debugging only)

ii              =   1;
hwait = waitbar(0,'Division in 3 parts...');
for i=1:length(permit_inds)
    img     =   squeeze(dataset(:,:,:,i)); 
    MSK     =   double(mask_fin(:,:,i));
    
    %% feature extraction
    img_hsv     =   rgb2hsv(img);
    img_cielab  =   applycform(img, reg); % eq. CIELAB
    
%     TLanti(i)   = uint16(fminbnd(@(x) sym_dissimilar_MSK(x,img_hsv,MSK,NBINs,30*SUBfac,0.5),search_range_H(1),search_range_H(2)));
%     BUsim(i)    = uint16(fminbnd(@(x) sym_similar_MSKLR(x,img_hsv(1:TLanti(i),:,:),MSK(1:TLanti(i),:),NBINs,20*SUBfac,0.7),search_range_W(1),search_range_W(2)));
%     LEGsim(i)   = uint16(fminbnd(@(x) sym_similar_MSKLR(x,img_hsv(TLanti(i)+1:end,:,:),MSK(TLanti(i)+1:end,:),NBINs,20*SUBfac,0.7),search_range_W(1),search_range_W(2)));
% 	if maskon
% 		HDanti(i)   = uint16(fminbnd(@(x) sym_dissimilar_MSKH(x,img_hsv,MSK,20),5,double(TLanti(i))));
% 	else
% 		HDanti(i)   = uint16(fminbnd(@(x) sym_dissimilar_MSKH2(x,BUsim(i),TLanti(i),img_hsv,MSK,NBINs,20*SUBfac,30*SUBfac),bord,double(TLanti(i))-bord));
% 	end

    TLanti(i)   = uint16(fminbnd(@(x) dissym_div(x,img_hsv,MSK,delta(1),alpha),search_range_H(1),search_range_H(2)));
    BUsim(i)    = uint16(fminbnd(@(x) sym_div(x,img_hsv(1:TLanti(i),:,:),MSK(1:TLanti(i),:),delta(2),alpha),search_range_W(1),search_range_W(2)));
    LEGsim(i)   = uint16(fminbnd(@(x) sym_div(x,img_hsv(TLanti(i)+1:end,:,:),MSK(TLanti(i)+1:end,:),delta(2),alpha),search_range_W(1),search_range_W(2)));
	HDanti(i)   = uint16(fminbnd(@(x) sym_dissimilar_MSKH(x,img_hsv,MSK,delta(1)),5,double(TLanti(i))));
   
%     TLanti(i)   = uint16(fminbnd(@(x) sym_dissimilar(x,img_cielab),search_range_H(1),search_range_H(2)));
%     BUsim(i)    = uint16(fminbnd(@(x) sym_similar(x,img_cielab(1:TLanti(i),:,:)),search_range_W(1),search_range_W(2)));
%     LEGsim(i)   = uint16(fminbnd(@(x) sym_similar(x,img_cielab(BUsim(i)+1:end,:,:)),search_range_W(1),search_range_W(2)));
%     HDanti(i)   = uint16(fminbnd(@(x) sym_dissimilar(x,img_cielab),1,double(TLanti(i))));

    if plotY
        subplot(maxplot/10,maxplot/4,ii);
        imagesc(img), axis image,axis off,hold on;
        plot([1,W],[TLanti(i),TLanti(i)],'r','LineWidth',3);
        plot([1,W],[HDanti(i),HDanti(i)],'r','LineWidth',3);
        plot([BUsim(i),BUsim(i)],[HDanti(i),TLanti(i)],'r','LineWidth',3);
        plot([LEGsim(i),LEGsim(i)],[TLanti(i)+1,H],'r','LineWidth',3);hold off;
        ii = ii+1;
    end
    if ii>maxplot && plotY
        ii = 1;
        pause,clf(gcf)
    end
	waitbar(i/length(permit_inds),hwait)
end
close(hwait)

%% Head detection
% if dethead
%     try
% 		load(['MAT/headdet_' dataname '_f' num2str(SUBfac) '_Exp' num2str(expnum) '.mat'])% head detection loading (det_final)
%         fprintf('And detected head LOADING... ')
%     catch
%         fprintf('And detected head COMPUTATION... ')
%         MeerHeadDetection;
%         save(['MAT/headdet_' dataname '_f' num2str(SUBfac) '_Exp' num2str(expnum) '.mat'],'det_final');     %% saving
%     end
% else % no head detection alg.
%     det_final = NaN*ones(length(permit_inds),4);
% end
det_final = NaN*ones(length(permit_inds),4); % no head det

%% Kernel-map computation
ii = 1;
hwait = waitbar(0,'Kernel map...');
for i=1:length(permit_inds)
    img     =   squeeze(dataset(:,:,:,i)); 
    
    img_hsv     =   rgb2hsv(img);
    tmp         =   img_hsv(:,:,3);
    tmp         =   histeq(tmp); % Color Equalization
    img_hsv     =   cat(3,img_hsv(:,:,1),img_hsv(:,:,2),tmp); % eq. HSV


    if ~any(isnan(det_final(i,:))) % NaN = head not found
        HEAD = img_hsv(1:HDanti(i),:,:);
        cntr = [det_final(i,1)+det_final(i,3)/2,det_final(i,2)+det_final(i,4)/2];
        HEADW = radial_gau_kernel(cntr,DIMW*3,size(HEAD,1),W);
    else
        HEADW = zeros(HDanti(i),W);
    end

    if (HDanti(i)+1 >= TLanti(i))
        HDanti(i) = HDanti(i) - 2;
    end
    
    UP = img_hsv(HDanti(i)+1:TLanti(i),:,:);
    UPW = gau_kernel(BUsim(i),varW,size(UP,1),W);

    DOWN = img_hsv(TLanti(i)+1:end,:,:);
    DOWNW = gau_kernel(LEGsim(i),varW,size(DOWN,1),W);

    MAP_KRNL{i} = [HEADW/max(HEADW(:));UPW/max(UPW(:));DOWNW/max(DOWNW(:))];
    if (H-size(MAP_KRNL{i})>=0)
        MAP_KRNL{i} = padarray(MAP_KRNL{i},H-size(MAP_KRNL{i},1),'replicate','post');
    else
        MAP_KRNL{i} = MAP_KRNL{i}(1:H,:);
    end

    if plotY
        subplot(maxplot*2/12,12,ii), imagesc(img), axis image,axis off,hold on;
        subplot(maxplot*2/12,12,ii+1),imagesc(MAP_KRNL{i}), axis image
    end


    if ~any(isnan(det_final(i,:)))
        HEADW = HEADW(:);
        head_det_flag(i) = 1;
    else
        head_det_flag(i) = 0;
    end
    
    ii=ii+2;
    clear HEADW DOWNW UPW

    if ii>maxplot && plotY
        ii = 1;
        pause,clf(h1)
	end
	waitbar(i/length(permit_inds),hwait)
end
head_det = det_final; % NaN = head not found
close(hwait)