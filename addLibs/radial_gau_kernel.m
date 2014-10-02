function map = radial_gau_kernel(c,varW,H,W)

[xx,yy] = meshgrid(1:W,1:H);
xx = xx(:);
yy = yy(:);
map = mvnpdf([xx,yy],[c(1),c(2)],[varW,0;0,varW]);
map = reshape(map,[H,W]);