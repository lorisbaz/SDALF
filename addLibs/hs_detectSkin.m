function skin = hs_detectSkin(im, h_range, s_range)

%read in image & convert to HSV
im1 = double(im);
hsv_im1 = rgb2hsv(im1);

%pull out H & S data into new variables
H = hsv_im1(:,:,1);
S = hsv_im1(:,:,2);

%prevents errors
I = H+S;
I(find(I==0))=Inf;

h = H;
s = S;

%targets skin by only selecting values within the rectangle skin range
skin = ((s>s_range(1)) & (s<s_range(2)) &(h>h_range(1)) & (h<h_range(2)));