function visualize_blobs(img,mvec,pvec,num)

display = 1;
img=im2double(img);


[rows,cols,ndim]=size(img);

if display,
    figure(num);clf
    subplot(1,3,1);image((img));axis image
    title('Input image');
end

bkgr=[0 0 0]'; % Paint on black background
subplot(1,3,2)
bimg=draw_blobs(mvec,pvec,rows,cols,bkgr);
subplot(1,3,2);image((bimg));axis image

if display
      
    for i=1:size(mvec,2)
        bimg=draw_blobs(mvec(:,i),pvec(:,i),rows,cols,bkgr);
        subplot(1,3,3);image((bimg));axis image  
        title('MSCR output');
        pause
    end
end