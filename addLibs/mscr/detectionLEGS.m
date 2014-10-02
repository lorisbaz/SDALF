function [mvec, pvec, bimg] = detectionLEGS(img,region,num, display)

%display = 1;
%img=imread(img);
%img = illuminant_normalization(img);
img = uint8(img);
img=im2double(img); img = img(region,:,:);
[rows,cols,ndim]=size(img);

p.min_margin=0.003; %  0.0015;  % Set margin parameter
p.ainc = 1.05;
p.min_size = 40;
p. filter_size = 5;
p.verbosefl=0;  % Uncomment this line to suppress text output

if display,
    figure(num);clf
    subplot(1,2,1);image((img));axis image
    title('Input image');
end

[mvec,pvec, arate,elist2]=detect_mscr(img,p);

bkgr=[1 0 0]'; % Paint on black background
[mvec, pvec, arate] = eliminate_blobs_inthemask(mvec,pvec,arate ,img);
bimg=draw_blobs(mvec,pvec,rows,cols,bkgr);


if display
    subplot(1,2,2);image((bimg));axis image
    title('MSCR output');
end