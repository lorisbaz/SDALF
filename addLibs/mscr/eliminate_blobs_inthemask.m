
function [mvec, pvec, arate] =eliminate_blobs_inthemask(mvec,pvec, arate, img)

[dummy,Nbl]=size(mvec);
M=mvec(2:3,:);

elim_list = [];

for k=1:Nbl,
    
    c = img(round(M(2,k)), round(M(1,k)),:);
    if (c(1) == 0 && c(2) == 0 && c(3)==0)
        elim_list = [elim_list k];
    end
end

for i=1:size(pvec,2)
    if all(pvec(:,i)<0.05)
        elim_list = [elim_list i];
    end
end

keep_list = setdiff(1:Nbl, elim_list);

if isempty(keep_list)
    keep_list = 1:Nbl;
end

mvec = mvec(:,keep_list);
pvec = pvec(:,keep_list);
arate = arate(:,keep_list);

