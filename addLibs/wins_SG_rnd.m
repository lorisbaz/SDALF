function wins = wins_SG_rnd(R, N,fac,var,SG,perc,choise)
% wins = wins_SG_rnd(R, N,fac,var,SG,perc,choise)
%
% This function extract randomly a set of N overlapping windows (or patches)
% from an image or a region. In addiction, the segmentation can be inserted
% from input., and the algorithm cancel out the windows out of the
% segmentation.
%
% Input:
% - R: bounding box of the region of interest
% - N: numbser of patches (default 100)
% - fac: patch dimension (default [10,10])
% - var: variance of the patch dimension (default 4)
% - SG: segmentation
% - choise: type of sampling ('rand' default, or 'unif')
% - perc: cancel the border patches if the FG is below this %
%
% Output:
% - wins: a set of N bounding box (dim. Nx4)

% Copyright: Loris Bazzani   Date: 10/06/2009




if nargin < 7
    choise  = 'rand';
elseif nargin < 6
    perc    = 0.5;
elseif nargin < 5
    SG    = ones(R(4),R(3));
elseif nargin < 4
    var     = 4;
elseif nargin < 3
    fac     = [10,10];
elseif nargin < 2
    N       = 100;
end

W   = R(3);
H   = R(4);

if W < fac(1) || H < fac(2)
    wins = []; % controllo che la finestra possa contenere almeno una patch
else
    
    dimmean = [fac(1) fac(2)];
    sigma   = [var,0;0,var];
    wins    = zeros(N,4);
    
    wins = []; count = 0;
    while(size(wins,1) < N  && count < 5)
        % randomly generating centers and dimensions
        switch choise
            case 'rand'
                Px  = uint16(rand([N-size(wins,1),1])*(R(1)+R(3))+R(1));
                Py  = uint16(rand([N-size(wins,1),1])*(R(2)+R(4))+R(2));
            case 'unif'
                Px  = uint16(unifrnd(R(1),R(1)+R(3),[N-size(wins,1),1]));
                Py  = uint16(unifrnd(R(2),R(2)+R(4),[N-size(wins,1),1]));
        end
        dim = uint16(fliplr(mvnrnd(dimmean,sigma,N-size(wins,1))));
        
        % sampling from kernel map discrete probability
        cntrs = [Py,Px]; % [y,x]
        tmpind = cntrs(:,1)<=size(SG,1);
        cntrs = cntrs(tmpind,:);
        dim = dim(tmpind,:);
        tmpind = cntrs(:,2)<=size(SG,2);
        cntrs = cntrs(tmpind,:);
        dim = dim(tmpind,:);
        lin_cntrs   = sub2ind(size(SG),cntrs(:,1),cntrs(:,2));
        rand_sel    = rand(length(lin_cntrs),1);
        dim         = dim(SG(lin_cntrs)>=rand_sel,:); % using kernel for eliminating low probability windows
        cntrs       = cntrs(SG(lin_cntrs)>=rand_sel,:);
        wnow        = uint16([fliplr(cntrs),dim]);
        
        %     imagesc(SG),axis image
        %     hold on, scatter(Px,Py), hold off
        %     hold on, scatter(cntrs(:,2),cntrs(:,1),'r'), hold off
        
        
        % excedeed indexes
        wnow(wnow(:,1)<R(1) & wnow(:,3)~=0,1)      = R(1);
        wnow(wnow(:,2)<R(2) & wnow(:,4)~=0,2)      = R(2);
        wnow(wnow(:,1)+wnow(:,3)>R(1)+R(3)& wnow(:,3)~=0,3) = R(3)-wnow(wnow(:,1)+wnow(:,3)>R(1)+R(3)& wnow(:,3)~=0,1);
        wnow(wnow(:,2)+wnow(:,4)>R(2)+R(4)& wnow(:,4)~=0,4) = R(4)-wnow(wnow(:,2)+wnow(:,4)>R(2)+R(4)& wnow(:,4)~=0,2);
        wnow = wnow(wnow(:,3)>0,:);
        wnow = wnow(wnow(:,4)>0,:);
        
        % eliminating the too little patches (at the image border)
        tmpind = ~(wnow(:,3)<=dimmean(1)-sigma(1,1));
        wnow = wnow(tmpind,:);
        tmpind = ~(wnow(:,4)<=dimmean(2)-sigma(2,2));
        wnow = wnow(tmpind,:);
        
        if isempty(wnow)
            count = count + 1;
        end
        
        indtake = [];
        for w = 1:size(wnow,1)
            if sum(sum(SG(wnow(w,2):wnow(w,2)+wnow(w,4),wnow(w,1):wnow(w,1)+wnow(w,3))))>perc*double(wnow(w,4)*wnow(w,3))
                indtake = [indtake;w];
            end
        end
        
        wins    = [wins;wnow(indtake,:)];
        
        
        % zero dimensions patches
        zerodimwin = sum(wins(:,3)<=0 | wins(:,4)<=0);
    end
    
    if ~isempty(wins)
        wins    = wins(1:min(size(wins,1),N),:);
    end
    
end