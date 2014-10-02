function v = vtrans2(A,n)

%% M. Brand version


[r,c] = size(A);

B = permute(reshape(A,n,r/n,c),[1,3,2]);

v = reshape(B,n*c,r/n);

