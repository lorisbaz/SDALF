function patches = wind_extract_rnd(img,SG,N,dim,perc)
%  patches = wind_extract_rnd(img,SG,N,dim,perc)
%
% This function extract randomly a set of N overlapping windows (or patches)
% from an image or a region. In addiction, the segmentation can be inserted
% from input., and the algorithm cancel out the windows out of the
% segmentation.
%
% Input:
% - img: original image
% - SG: segmentation
% - N: numbser of patches (default 100)
% - dim: wins dimension [x,y] or [width,height]
% - perc: cancel the border patches if the FG is below this %
%
% Output:
% - patches: a set of patches (imgs)

% Copyright: Loris Bazzani   Date: 11/02/2009

R = [1,1,size(img,2),size(img,1)];
W   = R(3);
H   = R(4);
wins    = zeros(N,4);
	
wins = [];
while(size(wins,1) < N)   
    % randomly generating centers and dimensions	
	Px  = uint16(rand([N-size(wins,1),1])*(R(1)+R(3))+R(1));
	Py  = uint16(rand([N-size(wins,1),1])*(R(2)+R(4))+R(2));

    
    % sampling from kernel map discrete probability
    cntrs = [Py,Px]; % [y,x]
	cntrs = cntrs + fliplr(repmat(uint16(dim),size(cntrs,1),1))/2; % centres
    tmpind = cntrs(:,1)<=size(SG,1);
    cntrs = cntrs(tmpind,:);
    tmpind = cntrs(:,2)<=size(SG,2);
    cntrs = cntrs(tmpind,:);
    lin_cntrs   = sub2ind(size(SG),cntrs(:,1),cntrs(:,2));
    rand_sel    = rand(length(lin_cntrs),1);
    cntrs       = cntrs(SG(lin_cntrs)>=rand_sel,:);
    wnow        = uint16([fliplr(cntrs)-repmat(uint16(dim),size(cntrs,1),1)/2,repmat(dim,size(cntrs,1),1)]);
    
%     imagesc(SG),axis image
%     hold on, scatter(Px,Py), hold off
%     hold on, scatter(cntrs(:,2),cntrs(:,1),'r'), hold off


    % excedeed indexes
    wnow(wnow(:,1)<=R(1) & wnow(:,3)~=0,1)      = R(1);
    wnow(wnow(:,2)<=R(2) & wnow(:,4)~=0,2)      = R(2);
    wnow(wnow(:,1)+wnow(:,3)>=R(1)+R(3)& wnow(:,3)~=0,3) = R(3)-wnow(wnow(:,1)+wnow(:,3)>=R(1)+R(3)& wnow(:,3)~=0,1);
    wnow(wnow(:,2)+wnow(:,4)>=R(2)+R(4)& wnow(:,4)~=0,4) = R(4)-wnow(wnow(:,2)+wnow(:,4)>=R(2)+R(4)& wnow(:,4)~=0,2);
    wnow = wnow(wnow(:,3)>0,:);
    wnow = wnow(wnow(:,4)>0,:);
    
    indtake = [];
    for w = 1:size(wnow,1)
        if sum(sum(SG(wnow(w,2):wnow(w,2)+wnow(w,4)+1,wnow(w,1):wnow(w,1)+wnow(w,3)+1)))>perc*double(wnow(w,4)*wnow(w,3))
            indtake = [indtake;w];
        end        
    end

    wins    = [wins;wnow(indtake,:)];
       
    
    % zero dimensions patches
    zerodimwin = sum(wins(:,3)<=0 | wins(:,4)<=0);
end

wins    = wins(1:N,:);

%% patches extraction
patches = cell(N,1);
for n = 1:N
	patches{n} = img(wins(n,2):wins(n,2)+wins(n,4)+1,wins(n,1):wins(n,1)+wins(n,3)+1,:);
end

