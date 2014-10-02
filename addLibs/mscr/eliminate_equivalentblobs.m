function [mn, pn] = eliminate_equivalentblobs(mvec, pvec)

eliminate = [];

for i=1:size(mvec,2)
    c1 = mvec(2:3,i);
    for j=i+1:size(mvec,2)

        c2 = mvec(2:3,j);
        if sqrt( sum((c1-c2).^2) ) < 10 && mvec(1,i)/mvec(1,j)> 0.6 ... 
                && mvec(1,i)/mvec(1,j) < 1.4 && sqrt( sum((pvec(:,i)-pvec(:,j)).^2) )< 0.1
            eliminate = [eliminate j];
        end
        
    end
    
end

tokeep = setdiff(1:size(mvec,2), eliminate);
mn = mvec(:,tokeep);
pn = pvec(:,tokeep);
