%% part-based weighted HSV histogram
hwait = waitbar(0,'Extract wHSV..');
for i=1:length(permit_inds)

    img         =   squeeze(dataset(:,:,:,i));
    img_hsv     =   rgb2hsv(img);
    tmp         =   img_hsv(:,:,3);
    tmp         =   histeq(tmp); % Color Equalization
    img_hsv     =   cat(3,img_hsv(:,:,1),img_hsv(:,:,2),tmp); % eq. HSV

    Map = MAP_KRNL{i};

    if ~any(isnan(head_det(i,:)))
        HEAD = img_hsv(1:HDanti(i),:,:);
        cntr = [head_det(i,1)+head_det(i,3)/2,head_det(i,2)+head_det(i,4)/2];
        HEADW = Map(1:HDanti(i),:);
        head_det_flag(i) = 1;
    else
        head_det_flag(i) = 0;
    end

    UP = img_hsv(HDanti(i)+1:TLanti(i),:,:);
    UPW = Map(HDanti(i)+1:TLanti(i),:);

    DOWN = img_hsv(TLanti(i)+1:end,:,:);
    DOWNW = Map(TLanti(i)+1:end,:);

    if head_det_flag(i)
        HEADW = HEADW(:);      
    end
    UPW = UPW(:);
    DOWNW = DOWNW(:);
    tmph0 = []; tmph2 = [];
    tmpup0 = []; tmpup2 = [];
    tmpdown0 = []; tmpdown2 = [];
    for ch =1:3
        if head_det_flag(i)
            rasterHEAD  =   HEAD(:,:,ch);
            tmph2     =   [tmph2; whistcY(rasterHEAD(:),HEADW,[0:1/(NBINs(ch)-1):1])];
        else
            tmph2     =   [tmph2; zeros(length([0:1/(NBINs(ch)-1):1]),1)];
        end

        rasterUP  =   UP(:,:,ch);
        tmpup2     =   [tmpup2; whistcY(rasterUP(:),UPW,[0:1/(NBINs(ch)-1):1])];

        rasterDOWN  =   DOWN(:,:,ch);
        tmpdown2     =   [tmpdown2; whistcY(rasterDOWN(:),DOWNW,[0:1/(NBINs(ch)-1):1])];

    end
    whisto2(:,i)=[tmph2',tmpup2',tmpdown2']';
    whisto2(find(isnan(whisto2))) = 0;

    clear HEADW DOWNW UPW
    waitbar(i/length(permit_inds),hwait);
end

close(hwait);