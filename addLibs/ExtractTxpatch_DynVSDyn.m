%% Maximally-CC textured patch computation
%  with clustering with Mean Shift on LBP histogram and patch position + 1NN
clear data_e max_txpatch
hwait = waitbar(0,'Extract Epitexture..');
ii = 1;
for i=1:length(dset)
		dd = uint16(length(dset(i).index)/2);
		rnddd = 1:length(dset(i).index); % NO random selection of the data 4 dynamic feature

		% --- dyn 1 --- %
		for j = 1:dd % for each view/2
			img         =   squeeze(dataset(:,:,:,dset(i).index(rnddd(j))));
			img_hsv     =   rgb2hsv(img);
			tmp         =   img_hsv(:,:,3);
			tmp         =   histeq(tmp); % Color Equalization
			img_hsv     =   cat(3,img_hsv(:,:,1),img_hsv(:,:,2),tmp);
			
			img_gray	= rgb2gray(img);
			
			Map = MAP_KRNL{i}.*double(mask_fin(:,:,dset(i).index(rnddd(j))));
			
			div = uint16([1,BUsim(dset(i).index(rnddd(j)))+1,TLanti(dset(i).index(rnddd(j)))+1;
				BUsim(dset(i).index(rnddd(j))),TLanti(dset(i).index(rnddd(j))),H]);
			for part = 2:size(div,2) % NO head
				
				img_part = img_gray(div(1,part):div(2,part),:);
				R = [1,1,fliplr(size(img_part))-1]; % [x,y,dimx,dimy] = [c,r,dimc,dimr] = [w,h,dimw,dimh]
				patches =  wins_SG_rnd(R,N,fac(part,:),var,Map(div(1,part):div(2,part),:),double(mask_fin(:,:,dset(i).index(rnddd(j)))));
				patches(:,1:2) =  uint16(patches(:,1:2)) + repmat([0,div(1,part)-1],size(patches,1),1);
				
				histo_entropy = zeros(size(patches,1),1);
				w_ncc = zeros(size(patches,1),1);
				for p = 1:size(patches,1)
					ptcs_t = cell(NTrans,1);
					patch = img_gray(patches(p,2):patches(p,2)+patches(p,4),patches(p,1):patches(p,1)+patches(p,3));
					
					ptcs_t = ptcs_yTrans(patch,NTrans); % patch transformations
					
					patch_op  =  [floor(patches(p,1)-(DIM_OP(1)-patches(p,3))/2) floor(patches(p,2)-(DIM_OP(2)-patches(p,4))/2) DIM_OP(1) DIM_OP(2)];
					patch_op  =	 bbox_exceedes(double(patch_op),H,W);
					img_local =  img_gray(patch_op(2):patch_op(2)+patch_op(4),patch_op(1):patch_op(1)+patch_op(3));
					
					% RGB hist entropy of the patch computation
					histo_entropy(p)  = entropy(img(patches(p,2):patches(p,2)+patches(p,4),patches(p,1):patches(p,1)+patches(p,3),1)) ...
						+ entropy(img(patches(p,2):patches(p,2)+patches(p,4),patches(p,1):patches(p,1)+patches(p,3),2)) ...
						+ entropy(img(patches(p,2):patches(p,2)+patches(p,4),patches(p,1):patches(p,1)+patches(p,3),3));
					
					% Normalized Cross Correlation (NCC)
					if histo_entropy(p)>Thresh_entr
						if plotY
							rectangle('Position',patches(p,:),'EdgeColor',[1 0 0]);
						end
						map_ncc = abs(normxcorr2(patch,img_local)); % (original patch)
						map_ncc = map_ncc(ceil((size(patch,1))/2):ceil(size(img_local,1)+(size(patch,1))/2)-1,...
							ceil((size(patch,2))/2):ceil(size(img_local,2)+(size(patch,2))/2)-1);
					else
						map_ncc = zeros(size(img_local));
					end
					for t=1:length(ptcs_t);
						if sum(ptcs_t{t}(:)>0) && histo_entropy(p)>Thresh_entr ...
								&& ~all(ptcs_t{t}(:) == ptcs_t{t}(1,1))
							
							tmp = abs(normxcorr2(ptcs_t{t},img_local)); % (transformed patch)
							tmp = tmp(ceil((size(ptcs_t{t},1))/2):ceil(size(img_local,1)+(size(ptcs_t{t},1))/2)-1,...
								ceil((size(ptcs_t{t},2))/2):ceil(size(img_local,2)+(size(ptcs_t{t},2))/2)-1);
							map_ncc = map_ncc + tmp;
							
						end
					end
					map_ncc = map_ncc/(NTrans+1); % normalization
					MNCC{p} = map_ncc;
					w_ncc(p) = sum(abs(map_ncc(:)))/numel(img_local);
					
					% CC thresholding and connected component analysis
					thr_map_cc = (map_ncc>Thresh_CC);
					area_ncctmp = regionprops(bwlabel(thr_map_cc),'Area');
					if ~isempty(area_ncctmp)
						area_ncc(p) = 1-max(struct2array(area_ncctmp))/numel(map_ncc); % max normalized area
					else
						area_ncc(p) = 0;
					end
					
					tmpPATCH{p} = img(patches(p,2):patches(p,2)+patches(p,4),patches(p,1):patches(p,1)+patches(p,3),:);
					
				end
				
				%% maximally-texturized patch computation (ranking by NCC of the first K high-entropy patches)
				[v,order]	= sort(histo_entropy,'descend'); % entropy sorting
				[v,order2]	= sort(area_ncc(order),'descend'); % NCC sorting
				ep = 1;
				for p = 1:size(patches,1)
					if histo_entropy(order(order2(p)))>Thresh_entr  ...
							&& w_ncc(order(order2(p)))>0
						epitome(j,part).patch{ep} = tmpPATCH{order(order2(p))};
						epitome(j,part).entr{ep}  = histo_entropy(order(order2(p)));
						epitome(j,part).w_ncc{ep} = area_ncc(order(order2(p)));
						epitome(j,part).pos{ep}   = patches(order(order2(p)),:);
						epitome(j,part).thrCC{ep} = MNCC{order(order2(p))};
						
						% LBP histogram
						lbph = lbp(epitome(j,part).patch{ep}, 1,8,mapping, 'h');
						
						% patch centroid
						ctr = double([epitome(j,part).pos{ep}(1)+epitome(j,part).pos{ep}(3)/2,...
							epitome(j,part).pos{ep}(2)+epitome(j,part).pos{ep}(4)/2]);
						
						epitome(j,part).data_e(:,ep) = [ctr';lbph']; % justapposition of data_e
						
						ep = ep + 1;
					end
				end
				try
					epitome(j,part);
				catch % if no epitome is found
					epitome(j,part).patch = [];
				end
				
				
			end
		end
		
		%% MS Clustering
		for part = 2:size(div,2) % NO head
			data_e			= [];
			inddata_e		= [];
			subinddata_e	= [];
			for j = 1:dd
				data_e		= [data_e,epitome(j,part).data_e];
				inddata_e	= [inddata_e,repmat(j,1,size(epitome(j,part).data_e,2))];
				subinddata_e= [subinddata_e,1:size(epitome(j,part).data_e,2)];
			end

			if ~isempty(inddata_e)
				[clust_res,centres]= MS_Clustering2(data_e',MSpar,'Uniform',[2,size(data_e,1)-2],2);
				K = max(clust_res);
				
				% 1-NN to select representative element
				if K~=0
					for k = 1:K
						inds = find(clust_res == k);
						if length(inds)>1
							D = squareform(pdist([data_e(:,inds)';centres(k,:)]));
							D = D(1:size(data_e(:,inds),2),size(data_e(:,inds),2)+1:end);
							[v minind] = min(D);
							rap(k) = inds(minind);
							dim_cluster(k) = length(inds);
						else
							rap(k) = -1; % eliminate single-element clusters
							dim_cluster(k) = -1;
						end
					end
				else
					rap = repmat(-1,1,size(data_e,2));
				end
			end
			
			ep = 1;
			try
				for k = 1:length(rap)
					if rap(k)>0
						max_txpatch(ii,part).patch{ep}	= epitome(inddata_e(rap(k)),part).patch{subinddata_e(rap(k))};
						max_txpatch(ii,part).entr{ep}	= epitome(inddata_e(rap(k)),part).entr{subinddata_e(rap(k))};
						max_txpatch(ii,part).w_ncc{ep}	= epitome(inddata_e(rap(k)),part).w_ncc{subinddata_e(rap(k))};
						max_txpatch(ii,part).pos{ep}	= epitome(inddata_e(rap(k)),part).pos{subinddata_e(rap(k))};
						max_txpatch(ii,part).lbph{ep}	= epitome(inddata_e(rap(k)),part).data_e(3:end,subinddata_e(rap(k))); % only lpb hist.
						max_txpatch(ii,part).numel{ep}	= dim_cluster(k);
						ep = ep + 1;
					end
				end
			catch % if no epitome is found
				max_txpatch(ii,part).patch	= [];
				max_txpatch(ii,part).entr	= [];
				max_txpatch(ii,part).w_ncc	= [];
				max_txpatch(ii,part).pos	= [];
				max_txpatch(ii,part).lbph	= [];
				max_txpatch(ii,part).numel	= [];
			end
			clear rap dim_cluster
		end
		ii = ii+1;
		clear epitome
		
		
		
		
		% --- dyn 2 --- %
		for j = dd+1:length(dset(i).index) % for each view/2
			img         =   squeeze(dataset(:,:,:,dset(i).index(rnddd(j))));
			img_hsv     =   rgb2hsv(img);
			tmp         =   img_hsv(:,:,3);
			tmp         =   histeq(tmp); % Color Equalization
			img_hsv     =   cat(3,img_hsv(:,:,1),img_hsv(:,:,2),tmp);
			
			img_gray	= rgb2gray(img);
			
			Map = MAP_KRNL{i}.*double(mask_fin(:,:,dset(i).index(rnddd(j))));
			
			div = uint16([1,BUsim(dset(i).index(rnddd(j)))+1,TLanti(dset(i).index(rnddd(j)))+1;
				BUsim(dset(i).index(rnddd(j))),TLanti(dset(i).index(rnddd(j))),H]);
			for part = 2:size(div,2) % NO head
				
				img_part = img_gray(div(1,part):div(2,part),:);
				R = [1,1,fliplr(size(img_part))-1]; % [x,y,dimx,dimy] = [c,r,dimc,dimr] = [w,h,dimw,dimh]
				patches =  wins_SG_rnd(R,N,fac(part,:),var,Map(div(1,part):div(2,part),:),double(mask_fin(:,:,dset(i).index(rnddd(j)))));
				patches(:,1:2) =  uint16(patches(:,1:2)) + repmat([0,div(1,part)-1],size(patches,1),1);
				
				histo_entropy = zeros(size(patches,1),1);
				w_ncc = zeros(size(patches,1),1);
				for p = 1:size(patches,1)
					ptcs_t = cell(NTrans,1);
					patch = img_gray(patches(p,2):patches(p,2)+patches(p,4),patches(p,1):patches(p,1)+patches(p,3));
					
					ptcs_t = ptcs_yTrans(patch,NTrans); % patch transformations
					
					patch_op  =  [floor(patches(p,1)-(DIM_OP(1)-patches(p,3))/2) floor(patches(p,2)-(DIM_OP(2)-patches(p,4))/2) DIM_OP(1) DIM_OP(2)];
					patch_op  =	 bbox_exceedes(double(patch_op),H,W);
					img_local =  img_gray(patch_op(2):patch_op(2)+patch_op(4),patch_op(1):patch_op(1)+patch_op(3));
					
					% RGB hist entropy of the patch computation
					histo_entropy(p)  = entropy(img(patches(p,2):patches(p,2)+patches(p,4),patches(p,1):patches(p,1)+patches(p,3),1)) ...
						+ entropy(img(patches(p,2):patches(p,2)+patches(p,4),patches(p,1):patches(p,1)+patches(p,3),2)) ...
						+ entropy(img(patches(p,2):patches(p,2)+patches(p,4),patches(p,1):patches(p,1)+patches(p,3),3));
					
					% Normalized Cross Correlation (NCC)
					if histo_entropy(p)>Thresh_entr
						if plotY
							rectangle('Position',patches(p,:),'EdgeColor',[1 0 0]);
						end
						map_ncc = abs(normxcorr2(patch,img_local)); % (original patch)
						map_ncc = map_ncc(ceil((size(patch,1))/2):ceil(size(img_local,1)+(size(patch,1))/2)-1,...
							ceil((size(patch,2))/2):ceil(size(img_local,2)+(size(patch,2))/2)-1);
					else
						map_ncc = zeros(size(img_local));
					end
					for t=1:length(ptcs_t);
						if sum(ptcs_t{t}(:)>0) && histo_entropy(p)>Thresh_entr ...
								&& ~all(ptcs_t{t}(:) == ptcs_t{t}(1,1))
							
							tmp = abs(normxcorr2(ptcs_t{t},img_local)); % (transformed patch)
							tmp = tmp(ceil((size(ptcs_t{t},1))/2):ceil(size(img_local,1)+(size(ptcs_t{t},1))/2)-1,...
								ceil((size(ptcs_t{t},2))/2):ceil(size(img_local,2)+(size(ptcs_t{t},2))/2)-1);
							map_ncc = map_ncc + tmp;
							
						end
					end
					map_ncc = map_ncc/(NTrans+1); % normalization
					MNCC{p} = map_ncc;
					w_ncc(p) = sum(abs(map_ncc(:)))/numel(img_local);
					
					% CC thresholding and connected component analysis
					thr_map_cc = (map_ncc>Thresh_CC);
					area_ncctmp = regionprops(bwlabel(thr_map_cc),'Area');
					if ~isempty(area_ncctmp)
						area_ncc(p) = 1-max(struct2array(area_ncctmp))/numel(map_ncc); % max normalized area
					else
						area_ncc(p) = 0;
					end
					
					tmpPATCH{p} = img(patches(p,2):patches(p,2)+patches(p,4),patches(p,1):patches(p,1)+patches(p,3),:);
					
				end
				
				%% maximally-texturized patch computation (ranking by NCC of the first K high-entropy patches)
				[v,order]	= sort(histo_entropy,'descend'); % entropy sorting
				[v,order2]	= sort(area_ncc(order),'descend'); % NCC sorting
				ep = 1;
				for p = 1:size(patches,1)
					if histo_entropy(order(order2(p)))>Thresh_entr  ...
							&& w_ncc(order(order2(p)))>0
						epitome(j,part).patch{ep} = tmpPATCH{order(order2(p))};
						epitome(j,part).entr{ep}  = histo_entropy(order(order2(p)));
						epitome(j,part).w_ncc{ep} = area_ncc(order(order2(p)));
						epitome(j,part).pos{ep}   = patches(order(order2(p)),:);
						epitome(j,part).thrCC{ep} = MNCC{order(order2(p))};
						
						% LBP histogram
						lbph = lbp(epitome(j,part).patch{ep}, 1,8,mapping, 'h');
						
						% patch centroid
						ctr = double([epitome(j,part).pos{ep}(1)+epitome(j,part).pos{ep}(3)/2,...
							epitome(j,part).pos{ep}(2)+epitome(j,part).pos{ep}(4)/2]);
						
						epitome(j,part).data_e(:,ep) = [ctr';lbph']; % justapposition of data_e
						
						ep = ep + 1;
					end
				end
				try
					epitome(j,part);
				catch % if no epitome is found
					epitome(j,part).patch = [];
				end
				
				
			end
		end
		
		%% MS Clustering
		for part = 2:size(div,2) % NO head
			data_e			= [];
			inddata_e		= [];
			subinddata_e	= [];
			for j = 1:dd
				data_e		= [data_e,epitome(j,part).data_e];
				inddata_e	= [inddata_e,repmat(j,1,size(epitome(j,part).data_e,2))];
				subinddata_e= [subinddata_e,1:size(epitome(j,part).data_e,2)];
			end

			if ~isempty(inddata_e)
				[clust_res,centres]= MS_Clustering2(data_e',MSpar,'Uniform',[2,size(data_e,1)-2],2);
				K = max(clust_res);
				
				% 1-NN to select representative element
				if K~=0
					for k = 1:K
						inds = find(clust_res == k);
						if length(inds)>1
							D = squareform(pdist([data_e(:,inds)';centres(k,:)]));
							D = D(1:size(data_e(:,inds),2),size(data_e(:,inds),2)+1:end);
							[v minind] = min(D);
							rap(k) = inds(minind);
							dim_cluster(k) = length(inds);
						else
							rap(k) = -1; % eliminate single-element clusters
							dim_cluster(k) = -1;
						end
					end
				else
					rap = repmat(-1,1,size(data_e,2));
				end
			end
			
			ep = 1;
			try
				for k = 1:length(rap)
					if rap(k)>0
						max_txpatch(ii,part).patch{ep}	= epitome(inddata_e(rap(k)),part).patch{subinddata_e(rap(k))};
						max_txpatch(ii,part).entr{ep}	= epitome(inddata_e(rap(k)),part).entr{subinddata_e(rap(k))};
						max_txpatch(ii,part).w_ncc{ep}	= epitome(inddata_e(rap(k)),part).w_ncc{subinddata_e(rap(k))};
						max_txpatch(ii,part).pos{ep}	= epitome(inddata_e(rap(k)),part).pos{subinddata_e(rap(k))};
						max_txpatch(ii,part).lbph{ep}	= epitome(inddata_e(rap(k)),part).data_e(3:end,subinddata_e(rap(k))); % only lpb hist.
						max_txpatch(ii,part).numel{ep}	= dim_cluster(k);
						ep = ep + 1;
					end
				end
			catch % if no epitome is found
				max_txpatch(ii,part).patch	= [];
				max_txpatch(ii,part).entr	= [];
				max_txpatch(ii,part).w_ncc	= [];
				max_txpatch(ii,part).pos	= [];
				max_txpatch(ii,part).lbph	= [];
				max_txpatch(ii,part).numel	= [];
			end
			clear rap dim_cluster
		end
		ii = ii+1;
		clear epitome
		
		waitbar(i/length(dset),hwait);
end
close(hwait);