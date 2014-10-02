%% Dynamic feature computation vs static
clear wHSVDyn
for i = 1:size(whisto2,2)
	wHSVDyn(i).hist = whisto2(:,i);
end
clear whisto2

%% wHSV Hist distances computation 
hh = waitbar(0,'Distances computation wHSV...');
for i=1:size(wHSVDyn,2)
    for j=1:size(wHSVDyn,2)
        if head_det_flag(i) ~= head_det_flag(j)
            whisttmp1 = wHSVDyn(i).hist(sum(NBINs):end); % cancel the head hist! % weight histogram
            whisttmp2 = wHSVDyn(j).hist(sum(NBINs):end);
        else
            whisttmp1 = wHSVDyn(i).hist;
            whisttmp2 = wHSVDyn(j).hist;
        end
        
        distance(i,j) = bhattacharyya(whisttmp1,whisttmp2); % weight histogram
    end
    waitbar(i*j/(size(wHSVDyn,2)^2),hh);
end
close(hh);
final_dist_hist = distance;    
