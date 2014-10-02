function bbox = bbox_exceedes(bbox,HIMG,WIMG)

top		= bbox(2);
left	= bbox(1);
bottom	= bbox(4) + top;
right	= bbox(3) + left;

top     = 1*(top<1) + HIMG*(top>HIMG) + top*(top>=1 && top<=HIMG);
left    = 1*(left<1) + WIMG*(left>WIMG) + left*(left>=1 && left<=WIMG);
bottom  = 1*(bottom<1) + HIMG*(bottom>HIMG) + bottom*(bottom>=1 && bottom<=HIMG);
right   = 1*(right<1) + WIMG*(right>WIMG) + right*(right>=1 && right<=WIMG);

bbox(2) = top;
bbox(1) = left;
bbox(4) = bottom - top;
bbox(3) = right - left;