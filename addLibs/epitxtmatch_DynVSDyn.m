%% epitext Hist distances computation 
hh = waitbar(0,'Distances computation Epitexture...');
distance = zeros(size(epitexture,2),size(epitexture,2));
for i=1:size(epitexture,2)
	hists_i = epitexture(i).hist;
    for j=1:size(epitexture,2)
		hists_j = epitexture(j).hist;

		dm = zeros(size(hists_i,2),1);
		for ii = 1:size(hists_i,2)
			d = zeros(size(hists_j,2),1);
			for jj = 1:size(hists_j,2)
				d(jj) = bhattacharyya(hists_i(:,ii),hists_j(:,jj));
			end
			dm(ii) = min(d); % min distances matching
		end		
		
		distance(i,j) = sum(dm(:)); % mean of the distances

    end
    waitbar(i*j/(size(epitexture,2)^2),hh);
end
close(hh);
dist_epitext = distance;   