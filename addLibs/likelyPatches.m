function ordPatches = likelyPatches( ep, pcl, patchSize, N )

epSize = [size(ep,1),size(ep,2)];
ep = cat(2,cat(1,ep,ep,ep),cat(1,ep,ep,ep),cat(1,ep,ep,ep));
ordPatches = zeros(patchSize(1),patchSize(2),3,N);

for n=1:N
    [i j] = ind2sub( size(pcl), find( pcl == max(pcl(:)) == 1) );
    ordPatches(:,:,:,n) = ep(epSize(1)+i+1:i+epSize(1)+patchSize(1),epSize(2)+j+1:j+epSize(2)+patchSize(2),:);
%     subplot(1,N,n), imagesc( ordPatches(:,:,:,n) ), axis off
    pcl( pcl == max(pcl(:)) ) = -Inf;
end