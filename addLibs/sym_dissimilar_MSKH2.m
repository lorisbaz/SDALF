function d = sym_dissimilar_MSKH2(x,bust,low,img,MSK,quanti_hist,loc1,loc2)

warning off;

x = uint16(x);
[H,W,chs] = size(img);
imgUP   = double(img(1:x,:,:));
imgDOWN = double(img(x:low,:,:));
MSK_U       = MSK(1:x,:);
MSK_D       = MSK(x:end,:);

window = [max(bust-loc1/2,1),max(x-loc2/2,1),min(bust+loc1/2,W),min(x+loc2/2,low)]; %[xmin,ymin,xmax,ymax]

%% HSV histogram
histUP =[];
for ch =1:3
    rst      =   imgUP(window(2):end,window(1):window(3),ch)*quanti_hist(ch);
    histUP      =   [histUP; whistcY(rst(:),raster(MSK_U(window(2):end,window(1):window(3))),1:quanti_hist(ch))];
end

histDOWN =[];
for ch =1:3
    rst      =   imgDOWN(1:window(4)-x,window(1):window(3),ch)*quanti_hist(ch);
    histDOWN    =   [histDOWN; whistcY(rst(:),raster(MSK_D(1:window(4)-x,window(1):window(3))),1:quanti_hist(ch))];
end

 d = abs(1-bhattacharyya(histUP,histDOWN));