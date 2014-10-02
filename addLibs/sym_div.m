function d = sym_div(x,img,MSK,loc,alpha)

warning off;

x = uint16(x);
[H,W,chs] = size(img);
imgL    = double(img(:,1:x,:));
imgR    = double(img(:,x:end,:));
MSK_L       = MSK(:,1:x);
MSK_R       = MSK(:,x:end);

localderL  =  max(x-loc,1):x;
localderR  =  1:min(loc+1,size(MSK_R,2));
[dunny indM] = min([length(localderL),length(localderR)]); % different dim,
if indM == 1
	dimLoc = length(localderL);
elseif indM == 2
	dimLoc = length(localderR);
end
indexes = 1:dimLoc;


if dimLoc ~= 0
	% local imgs
	imgLloc(:,:,1) = imgL(:,indexes,1);
	imgLloc(:,:,2) = imgL(:,indexes,2);
	imgLloc(:,:,3) = imgL(:,indexes,3);
	imgRloc(:,:,1) = fliplr(imgR(:,indexes,1));% img flipping
	imgRloc(:,:,2) = fliplr(imgR(:,indexes,2));
	imgRloc(:,:,3) = fliplr(imgR(:,indexes,3));
	
	% FG map 
% 	MSK_Lloc = MSK_L(:,indexes);
% 	MSK_Rloc = fliplr(MSK_R(:,indexes)); % flipping
% 	[i,j] = find(MSK_Lloc >= 1 & MSK_Rloc >= 1);
% 	ind = sub2ind(size(MSK_Lloc),i,j);
% 	
% 	if ~isempty(ind)	% distance	
		d = (alpha)*sqrt(sum((imgLloc(:)-imgRloc(:)).^2))/dimLoc + (1-alpha)*abs(sum(MSK_R(:))-sum(MSK_L(:)))/max([numel(MSK_R),numel(MSK_L)]);
% 	else
% 		d = 1;
% 	end
else
	d = 1;
end