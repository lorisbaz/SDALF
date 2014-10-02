function d = sym_dissimilar_MSK(x,img,MSK,quanti_hist,loc,alpha)

warning off;

x = uint16(x);
[H,W,chs] = size(img);
imgUP   = double(img(1:x,:,:));
imgDOWN = double(img(x:end,:,:));
MSK_U       = MSK(1:x,:);
MSK_D       = MSK(x:end,:);

localderU  =  max(x-loc,1):x;
localderD  =  1:min(loc+1,size(MSK_D,1));

%% HSV histogram
histUP =[];
for ch =1:3
    rst      =   imgUP(localderU,:,ch)*quanti_hist(ch);
    histUP      =   [histUP; whistcY(double(rst(:)),raster(MSK_U(localderU,:)),1:quanti_hist(ch))];
end

histDOWN =[];
for ch =1:3
    rst      =   imgDOWN(localderD,:,ch)*quanti_hist(ch);
    histDOWN    =   [histDOWN; whistcY(double(rst(:)),raster(MSK_D(localderD,:)),1:quanti_hist(ch))];
end

% d1 = norm_entropy(imgUP(logical(MSK_U)));
% d2 = norm_entropy(imgDOWN(logical(MSK_D)));
if sum(sum(MSK_D(localderD,:)))~=0 && sum(sum(MSK_U(localderU,:)))~=0
    d = (alpha)*abs(1-bhattacharyya(histUP,histDOWN)) + (1-alpha)*abs(sum(MSK_U(:))-sum(MSK_D(:)))/max([numel(MSK_U),numel(MSK_D)]);
    % dissimilarity in appearance + similarity in FG shape
else
     d = 1;
end
% d = -d;
%bhattacharyya(histUP,histDOWN); %+ 