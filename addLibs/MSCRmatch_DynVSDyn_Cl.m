%% Dynamic feature computation
switch dataname
	case {'ETHZ1','ETHZ2','ETHZ3'}
		clear BlobsDyn
		ii = 1;
		for i = 1:length(dset) % for each ped
			
			% dyn 1
			BlobsDyn(ii).mvec = [];
			BlobsDyn(ii).pvec = [];
			dd = uint16(length(ped(i).dynmodels)/2);
			rnddd = ped(i).rnddd; 
			for j = 1:dd % for each view/2
				ind = find(ped(i).dynmodels(rnddd(j)) == dset(i).globalindex); % we are searching for the local index (in permits_ind)
				BlobsDyn(ii).mvec = [BlobsDyn(ii).mvec,Blobs(dset(i).index(ind)).mvec];
				BlobsDyn(ii).pvec = [BlobsDyn(ii).pvec,Blobs(dset(i).index(ind)).pvec];
			end
			ii = ii+1;
			
			% dyn 2
			BlobsDyn(ii).mvec = [];
			BlobsDyn(ii).pvec = [];
			for j = dd+1:length(ped(i).dynmodels) % for each view/2
				ind = find(ped(i).dynmodels(rnddd(j)) == dset(i).globalindex); % we are searching for the local index (in permits_ind)
				BlobsDyn(ii).mvec = [BlobsDyn(ii).mvec,Blobs(dset(i).index(ind)).mvec];
				BlobsDyn(ii).pvec = [BlobsDyn(ii).pvec,Blobs(dset(i).index(ind)).pvec];
			end
			ii = ii+1;
			
		end
	case 'iLIDS'
		clear BlobsDyn
		ii = 1;
		for i = 1:length(dset) % for each ped
			
			% dyn 1
			BlobsDyn(ii).mvec = [];
			BlobsDyn(ii).pvec = [];
			dd = uint16(length(ped(i).dynmodels)/2);
			rnddd = ped(i).rnddd; 
			for j = 1:dd % for each view/2
				ind = find(ped(i).dynmodels(rnddd(j)) == dset(i).globalindex); % we are searching for the local index (in permits_ind)
				BlobsDyn(ii).mvec = [BlobsDyn(ii).mvec,Blobs(dset(i).index(ind)).mvec];
				BlobsDyn(ii).pvec = [BlobsDyn(ii).pvec,Blobs(dset(i).index(ind)).pvec];
			end
			ii = ii+1;
			
			% dyn 2
			BlobsDyn(ii).mvec = [];
			BlobsDyn(ii).pvec = [];
			for j = dd+1:length(ped(i).dynmodels) % for each view/2
				ind = find(ped(i).dynmodels(rnddd(j)) == dset(i).globalindex); % we are searching for the local index (in permits_ind)
				BlobsDyn(ii).mvec = [BlobsDyn(ii).mvec,Blobs(dset(i).index(ind)).mvec];
				BlobsDyn(ii).pvec = [BlobsDyn(ii).pvec,Blobs(dset(i).index(ind)).pvec];
			end
			ii = ii+1;	
		end
end

%% MSCR distances computation
hh = waitbar(0,'Image features computation...');
bestk   = zeros(length(BlobsDyn),1);
bestmu  = zeros(5,length(BlobsDyn));
bestcov = zeros(5,5,length(BlobsDyn));
for i=1:length(BlobsDyn)
    %%---- MSE analysis
    mser    =   BlobsDyn(i);

    [keep_list,BB]= blobs_inside_image(mser.mvec,[H,W]); % used to compute BB
    keep_list = [1:size(mser.mvec,2)];

	for j=1:length(keep_list);
		centroid(:,j) = BB(1:2,1,j)./[W;H];
		colour(:,j) = applycform(mser.pvec(:,j)', reg);
	end

	% blobs clustering  for dynamic feature!
	data = [centroid;colour];
	[cl(i).bestk,bestpp,cl(i).bestmu,cl(i).bestcov,dl,countf] = mixtures4(data,kmin,min(kmax,size(data,2)),regularize,th,covoption);
	num_data = size(data,2);
	cl(i).label	 = zeros(num_data,1);
	ll		 = zeros(num_data,1);
	for j = 1:num_data
		vote = zeros(cl(i).bestk,1);
		for k=1:cl(i).bestk
			vote(k)= mvnpdf(data(:,j),cl(i).bestmu(:,k),cl(i).bestcov(:,:,k));
		end
		[ll(j), cl(i).label(j)] = max(vote);
	end
	
    dataf(i).Mmvec=centroid;
    dataf(i).Mpvec=colour;

    waitbar(i/length(BlobsDyn),hh);
	
	clear centroid colour

end
close(hh);


%%
hh = waitbar(0,'Distances computation MSCR...');
for i = 1:length(dataf)
	
	Mmvec{1} = dataf(i).Mmvec;
	Mpvec{1} = dataf(i).Mpvec;
	
	num(1)=size(Mmvec{1},2);
	
	dist_y_n	= cell(length(permit_inds),1);
	dist_color_n= cell(length(permit_inds),1);
	
	for j = 1:length(dataf)
		
		if i == j
			
			final_dist_y(i,j)		= 0;
			final_dist_color(i,j)	= 0;
		else
			Mmvec{2} = dataf(j).Mmvec;
			Mpvec{2} = dataf(j).Mpvec;
			
			if mod(i,2) == 1 && mod(j,2) == 0 % dyn vs dyn
				num(2)= size(Mmvec{2},2); %length(Mindeces{2});
				[thrash,max_ind(i,j)]= max([num(1),num(2)]);
				max_i = mod(max_ind(i,j),2)+1;%max_ind;
				min_i = max_ind(i,j);%mod(max_ind,2)+1;
				max_info = num(max_i);
				min_info = num(min_i);
				
				dist_y = 1000*ones(min_info, max_info);
				dist_color = 1000*ones(min_info, max_info);
				
% 				% Cluster Association
% 				Dass = zeros(cl(i).bestk,cl(j).bestk);
% 				for cli = 1:cl(i).bestk
% 					for clj = 1:cl(j).bestk
% 						Dass(cli,clj) = mahalan(cl(i).bestmu(:,cli),cl(j).bestmu(:,clj),cl(j).bestcov(:,:,clj)) + ...
% 							mahalan(cl(j).bestmu(:,clj),cl(i).bestmu(:,cli),cl(i).bestcov(:,:,cli)); % d(mu2;mu1,sigma1) + d(mu1;mu2,sigma2)
% 					end
% 				end
% 				[trash max_i_b] = max([cl(i).bestk,cl(j).bestk]);
% 				[trash min_i_b] = min([cl(i).bestk,cl(j).bestk]);
% 				[trash indAss] = min(Dass,[],min_i_b); % minimize along the class with the biggest N of blobs 
% 				indAss = [[1:length(indAss)]',reshape(indAss,length(indAss),1)]; % correspondences matrix
% 				% find blobs ids for each correspondence found
% 				for as = 1:size(indAss,1)
% 					bl_into_cl1 = find(cl(i).label == indAss(as,max_i_b)); 
% 					bl_into_cl2 = find(cl(j).label == indAss(as,min_i_b));
% 					if min_i == 1 % swap indexes		
% 						for k=1:length(bl_into_cl2)
% 							for h=1:length(bl_into_cl1)
% 								dist_y(bl_into_cl1(h),bl_into_cl2(k)) = abs(Mmvec{max_i}(2,bl_into_cl2(k))-Mmvec{min_i}(2,bl_into_cl1(h)));
% 								dist_color(bl_into_cl1(h),bl_into_cl2(k)) = sqrt(sum((Mpvec{max_i}(:,bl_into_cl2(k))-Mpvec{min_i}(:,bl_into_cl1(h))).^2));
% 							end
% 						end
% 					else % min_i == 2
% 						for k=1:length(bl_into_cl2)
% 							for h=1:length(bl_into_cl1)
% 								dist_y(bl_into_cl2(k),bl_into_cl1(h)) = abs(Mmvec{max_i}(2,bl_into_cl1(h))-Mmvec{min_i}(2,bl_into_cl2(k)));
% 								dist_color(bl_into_cl2(k),bl_into_cl1(h)) = sqrt(sum((Mpvec{max_i}(:,bl_into_cl1(h))-Mpvec{min_i}(:,bl_into_cl2(k))).^2));
% 							end
% 						end
% 					end
% 				end
				
				
				for k=1:max_info % smallest
					% cluster selection
					vote = zeros(cl(i).bestk,1);
					for kl = 1:cl(i).bestk
						vote(kl)= mvnpdf([Mmvec{max_i}(:,k);Mpvec{max_i}(:,k)],cl(i).bestmu(:,kl),cl(i).bestcov(:,:,kl));
					end
					[trash NNind] = max(vote); % Nearest Neighbor data association
					bl_into_cl = find(cl(i).label == NNind); % blobs into selected cluster
					
					for h=1:length(bl_into_cl) % biggest
						
						dist_y(bl_into_cl(h),k) = abs(Mmvec{max_i}(2,k)-Mmvec{min_i}(2,bl_into_cl(h)));
						dist_color(bl_into_cl(h),k) = sqrt(sum((Mpvec{max_i}(:,k)-Mpvec{min_i}(:,bl_into_cl(h))).^2));
						
					end
				end
				
				
				% check outliers
				ref_y = min(dist_y); me_ref_y = mean(ref_y); std_ref_y = std(ref_y);
				ref_color = min(dist_color); me_ref_color = mean(ref_color); std_ref_color = std(ref_color);
				good = find((ref_y<(me_ref_y+3.5*std_ref_y))&(ref_color<(me_ref_color+3.5*std_ref_color)));
				
				max_useful_info = length(good);
				dist_y2 = dist_y(:,good);
				dist_color2 = dist_color(:,good);
				
				%%normalize
				dist_y_n{j} = dist_y./max(dist_y2(:));
				dist_color_n{j} = dist_color./max(dist_color2(:));
				
				%%distance computation
				totdist_n = ( pyy*dist_y_n{j}(:,good)  + pcc*dist_color_n{j}(:,good) ) ;
				
				%%Minimization
				[min_dist,matching] = min(totdist_n);
				
				useful_i = sub2ind(size(totdist_n),[matching]',[1:max_useful_info]');
				final_dist_y(i,j)  = sum(dist_y2(useful_i))/max_useful_info;
				final_dist_color(i,j) = sum(dist_color2(useful_i))/max_useful_info;
			else
				final_dist_y(i,j)		= 1;
				final_dist_color(i,j)	= 1;
			end
			
		end
	end
	waitbar(i*j/length(dataf)^2,hh);
	clear dist_y;
	clear dist_color;
	clear dist_color_n;
	clear dist_y_n;
	
end
close(hh);


for i=1:length(dataf)
	final_dist_y(i,:) =  final_dist_y(i,:)./max(final_dist_y(i,:));
	final_dist_color(i,:) = final_dist_color(i,:)./max(final_dist_color(i,:));
	final_mscr_dist(i,:) = (pyy*final_dist_y(i,:) +  pcc*final_dist_color(i,:) ); % +  final_dist_hist(i,:);
    
end
