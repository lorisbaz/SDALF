%% Dynamic feature computation
switch dataname
	case {'ETHZ1','ETHZ2','ETHZ3'}
		clear wHSVDyn
		ii = 1; ddtmp = 1;
		for i = 1:length(dset) % for each ped
			
			% dyn 1
			wHSVDyn(ii).hist = [];
			dd = uint16(length(ped(i).dynmodels)/2);
			rnddd = ped(i).rnddd; 
			for j = 1:dd % for each view/2
				ind = find(ped(i).dynmodels(rnddd(j)) == dset(i).globalindex); % we are searching for the local index (in permits_ind)
				wHSVDyn(ii).hist  =  [wHSVDyn(ii).hist,whisto2(:,dset(i).index(ind))];
			end
			ii = ii+1;
			
			% dyn 2
			wHSVDyn(ii).hist = [];
			for j = dd+1:length(ped(i).dynmodels) % for each view/2
				ind = find(ped(i).dynmodels(rnddd(j)) == dset(i).globalindex); % we are searching for the local index (in permits_ind)
				wHSVDyn(ii).hist  =  [wHSVDyn(ii).hist,whisto2(:,dset(i).index(ind))];
			end
			ii = ii+1;
			
		end
	case 'iLIDS'
		clear wHSVDyn
		ii = 1; ddtmp = 1;
		for i = 1:length(dset) % for each ped
			
			% dyn 1
			wHSVDyn(ii).hist = [];
			dd = uint16(length(ped(i).dynmodels)/2);
			rnddd = ped(i).rnddd; 
			for j = 1:dd % for each view/2
				ind = find(ped(i).dynmodels(rnddd(j)) == dset(i).globalindex); % we are searching for the local index (in permits_ind)
				wHSVDyn(ii).hist  =  [wHSVDyn(ii).hist,whisto2(:,dset(i).index(ind))];
			end
			ii = ii+1;
			
			% dyn 2
			wHSVDyn(ii).hist = [];
			for j = dd+1:length(ped(i).dynmodels) % for each view/2
				ind = find(ped(i).dynmodels(rnddd(j)) == dset(i).globalindex); % we are searching for the local index (in permits_ind)
				wHSVDyn(ii).hist  =  [wHSVDyn(ii).hist,whisto2(:,dset(i).index(ind))];
			end
			ii = ii+1;
			
		end
end
	
%% wHSV Hist distances computation 
hh = waitbar(0,'Distances computation wHSV...');
distance = zeros(size(wHSVDyn,2),size(wHSVDyn,2));
for i=1:size(wHSVDyn,2)
	hists_i = wHSVDyn(i).hist;
    for j=1:size(wHSVDyn,2)
		hists_j = wHSVDyn(j).hist;
        if head_det_flag(i) ~= head_det_flag(j)
            histtmp_i = hists_i(sum(NBINs):end,:); % cancel the head hist! % weight histogram
            histtmp_j = hists_j(sum(NBINs):end,:);
		else			
            histtmp_i = hists_i;
            histtmp_j = hists_j;	
		end

		dm = zeros(size(histtmp_i,2),1);
		for ii = 1:size(histtmp_i,2)
			d = zeros(size(histtmp_j,2),1);
			for jj = 1:size(histtmp_j,2)
				d(jj) = bhattacharyya(histtmp_i(:,ii),histtmp_j(:,jj));
			end
			dm(ii) = min(d); % min distances matching
		end		
		
		distance(i,j) = sum(dm(:)); % mean of the distances

    end
    waitbar(i*j/(size(wHSVDyn,2)^2),hh);
end
close(hh);
final_dist_hist = distance;   