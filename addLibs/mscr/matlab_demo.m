%
% Simple demo of Matlab MEX interface
%
%img=imread('puzzle.png');
img = imread('../../VIPeR/cam_a/000_45.bmp');
img=im2double(img);
[rows,cols,ndim]=size(img);


load ../Masks
B = M00(:,:,1);


p.min_margin=0.003; %  0.0015;  % Set margin parameter
p.ainc = 1.05;
p.min_size = 40;
%p.verbosefl=0;  % Uncomment this line to suppress text output

figure;clf
subplot(1,3,1);image(img);axis image
title('Input image');

subplot(1,3,2);imagesc(B);axis image
title('Mask image');

[mvec,pvec]=detect_mscr_masked(img,B,p);

bkgr=[0 0 1]'; % Paint on back background
bimg=draw_blobs(mvec,pvec,rows,cols,bkgr);
subplot(1,3,3);image(bimg);axis image
title('MSCR output');