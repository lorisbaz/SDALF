function d = sym_dissimilar_MSKH(x,img,MSK,loc)

warning off;

x = uint16(x);
[H,W,chs] = size(img);
imgUP   = double(img(1:x,:,:));
imgDOWN = double(img(x:end,:,:));
MSK_U       = MSK(1:x,:);
MSK_D       = MSK(x:end,:);

localderU  =  max(x-loc,1):x;
localderD  =  1:min(loc+1,size(MSK_D,1));

% if sum(sum(MSK_D(localderD,:)))~=0 || sum(sum(MSK_U(localderU,:)))~=0
    d = - abs(sum(sum(MSK_U(localderU,:)))-sum(sum(MSK_D(localderD,:))));
% else
%     d = 0;
% end
