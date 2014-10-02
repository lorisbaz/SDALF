function d = dissym_div_UP(x,img,MSK,loc,alpha)

warning off;

x = uint16(x);
[H,W,chs] = size(img);
imgUP   = double(img(1:x,:,:));
imgDOWN = double(img(x:end,:,:));
MSK_U       = MSK(1:x,:);
MSK_D       = MSK(x:end,:);

localderU  =  max(x-loc,1):x;
localderD  =  1:min(loc+1,size(MSK_D,1));
[dunny indM] = min([length(localderU),length(localderD)]); % different dim,
if indM == 1
	dimLoc = length(localderU);
elseif indM == 2
	dimLoc = length(localderD);
end
indexes = 1:dimLoc;


if dimLoc ~= 0
	% local imgs
	imgUPloc(:,:,1) = imgUP(indexes,:,1);
	imgUPloc(:,:,2) = imgUP(indexes,:,2);
	imgUPloc(:,:,3) = imgUP(indexes,:,3);
	imgDWloc(:,:,1) = fliplr(imgDOWN(indexes,:,1));% img flipping
	imgDWloc(:,:,2) = fliplr(imgDOWN(indexes,:,2));
	imgDWloc(:,:,3) = fliplr(imgDOWN(indexes,:,3));
	
	% FG map 
% 	MSK_Lloc = MSK_L(:,indexes);
% 	MSK_Rloc = fliplr(MSK_R(:,indexes)); % flipping
% 	[i,j] = find(MSK_Lloc >= 1 & MSK_Rloc >= 1);
% 	ind = sub2ind(size(MSK_Lloc),i,j);
% 	
% 	if ~isempty(ind)	% distance	
		d = (alpha)*(1-sqrt(sum((imgUPloc(:)-imgDWloc(:)).^2))/dimLoc) + (1-alpha)*(1-abs(sum(MSK_U(:))-sum(MSK_D(:)))/max([numel(MSK_U),numel(MSK_D)]));
% 	else
% 		d = 1;
% 	end
else
	d = 1;
end