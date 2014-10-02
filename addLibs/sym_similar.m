function d = sym_similar(x,img)

warning off;

[H,W,chs] = size(img);
imgLEFT   = img(:,1:x,:);
imgRIGHT  = img(:,x:end,:);

% if length(1:x)>length(x:W) % pixels cutting (l'alternativa e` il padding)
%     imgLEFT = imgLEFT(:,length(1:x)-length(x:W)+1:end,:);
% elseif length(1:x)<length(x:W)
%     imgRIGHT = imgRIGHT(:,1:length(1:x),:);
% end

if length(1:x)>length(x:W) %  padding
    imgRIGHT = padarray(imgRIGHT,[0,length(1:x)-length(x:W)],'replicate','post');
elseif length(1:x)<length(x:W)
    imgLEFT = padarray(imgLEFT,[0,length(x:W)-length(1:x)],'replicate','pre');
end

% img flipping
imgRIGHT(:,:,1) = fliplr(imgRIGHT(:,:,1));
imgRIGHT(:,:,2) = fliplr(imgRIGHT(:,:,2));
imgRIGHT(:,:,3) = fliplr(imgRIGHT(:,:,3));

% subplot(121),imagesc(imgLEFT),axis image
% subplot(122),imagesc(imgRIGHT),axis image

 
d = sqrt(sum((imgLEFT(:)-imgRIGHT(:)).^2)); % max the similarity dist.

% quanti_hist = 16;
% tmp =[];
% for ch =1:chs
%     raster  =   imgLEFT(:,:,ch);
%     tmp     =   [tmp, hist(raster(:),quanti_hist)];
% end
% histUP = tmp;
% 
% tmp =[];
% for ch =1:chs
%     raster  =   imgRIGHT(:,:,ch);
%     tmp     =   [tmp, hist(raster(:),quanti_hist)];
% end
% histDOWN = tmp;
% 
% d = bhattacharyya(histUP,histDOWN);