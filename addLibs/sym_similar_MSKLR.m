function d = sym_similar_MSKLR(x,img,MSK,quanti_hist,loc,alpha)

warning off;

x = uint16(x);
[H,W,chs] = size(img);
imgL    = double(img(:,1:x,:));
imgR    = double(img(:,x:end,:));
MSK_L       = MSK(:,1:x);
MSK_R       = MSK(:,x:end);

localderL  =  max(x-loc,1):x;
localderR  =  1:min(loc+1,size(MSK_R,2));

%% HSV histogram
histL =[];
for ch =1:3
    rst      =   imgL(:,localderL,ch)*quanti_hist(ch);
    histL    =   [histL; whistcY(rst(:),raster(MSK_L(:,localderL)),1:quanti_hist(ch))];
end

histR =[];
for ch =1:3
    rst      =   imgR(:,localderR,ch)*quanti_hist(ch);
    histR    =   [histR; whistcY(rst(:),raster(MSK_R(:,localderR)),1:quanti_hist(ch))];
end

% if sum(sum(MSK_L(:,localderR)))~=0 && sum(sum(MSK_R(:,localderL)))~=0
d = (alpha)*bhattacharyya(histL,histR) + (1-alpha)*abs(sum(MSK_R(:))-sum(MSK_L(:)))/max([numel(MSK_R),numel(MSK_L)]);
    % similarity in appearance + simmmetry of FG shape
% else
%      d = 1;
% end
%