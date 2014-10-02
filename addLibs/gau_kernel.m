function map = gau_kernel(x,varW,H,W)

map = repmat(normpdf([1:W],double(x),varW),[H,1]);