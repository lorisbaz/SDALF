function [mvec, pvec, bimg] = detection(img, mask, region,p,num, display)

img=im2double(img); img = img(region,:,:);
mask = mask(region,:);

[rows,cols,ndim]=size(img);


if display,
    figure(num);clf
    subplot(1,2,1);image((img));axis image
    title('Input image');
end

[mvec,pvec, arate,elist2]=detect_mscr_masked(img,mask,p);
pvec = pvec/256;

bkgr=[1 0 0]'; % Paint on black background
[mvec, pvec] = eliminate_equivalentblobs(mvec, pvec);
bimg=draw_blobs(mvec,pvec,rows,cols,bkgr);


if display
      subplot(1,2,2);image((bimg));axis image
%     for i=1:size(mvec,2)
%         bimg=draw_blobs(mvec(:,i),pvec(:,i),rows,cols,bkgr);
%         subplot(1,2,2);image((bimg));axis image  
%         title('MSCR output');
%         pause
%     end
end