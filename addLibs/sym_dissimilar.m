function d = sym_dissimilar(x,img)

warning off;

[H,W,chs] = size(img);
imgUP   = img(1:x,:,:);
imgDOWN = img(x:end,:,:);

if length(1:x)>length(x:H) % pixels cutting (l'alternativa e` il padding)
    imgUP = imgUP(length(1:x)-length(x:H)+1:end,:,:);
elseif length(1:x)<length(x:H)
    imgDOWN = imgDOWN(1:length(1:x),:,:);
end

% if length(1:x)>length(x:H) %  padding
%     imgDOWN = padarray(imgDOWN,[length(1:x)-length(x:H),0],'post');
% elseif length(1:x)<length(x:H)
%     imgUP = padarray(imgUP,[length(x:H)-length(1:x),0],'pre');
% end

% subplot(211),imagesc(imgUP),axis image
% subplot(212),imagesc(imgDOWN),axis image
d = -sqrt(sum((imgUP(:)-imgDOWN(:)).^2)); % min the similarity dist.

% quanti_hist = 16;
% tmp =[];
% for ch =1:chs
%     raster  =   imgUP(:,:,ch);
%     tmp     =   [tmp, hist(raster(:),quanti_hist)];
% end
% histUP = tmp;
% 
% tmp =[];
% for ch =1:chs
%     raster  =   imgUP(:,:,ch);
%     tmp     =   [tmp, hist(raster(:),quanti_hist)];
% end
% histDOWN = tmp;
% 
% d = bhattacharyya(histUP,histDOWN);