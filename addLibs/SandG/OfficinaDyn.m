%% Dynamic model test 0
clear all, clc
load TMP_gauss_clust

avoid_inds  = [11,14,33,42,60,67,68,91,92,96,127,128,147,148,149,165,195,210,...
	212,217,218,236,248,259,260,290,303,310,328,427,439,39,110,120,166,...
	167,178,338];

dir_e   =   '/Users/bazzani/Pictures/i-LIDS_Pedestrian'; % '/Volumes/Untitled 1/datasets/Tagging/i-LIDS_Pedestrian';

%% Init
dir_e_list = dir(strcat(dir_e,'/*.jpg'));
n_img       =   size(dir_e_list,1);
permit_inds = setdiff(1:n_img,avoid_inds);
SUBfac  = 1;    % subsampling factor (for different-resolution experiments)
H = 128*SUBfac; W = 64*SUBfac; % NORMALIZED dimensions
reg  	=   makecform('srgb2lab');
pyy = 0.4;
pcc = 0.6;
pee = 0.7;


kmin=1;
kmax=5;
regularize = 0;
th = 1e-4;
covoption=1;


addpath('../addLibs/mscr/');
load TMP_gauss_clust
    
SetDataset;   
	
%% Dynamic feature building
ii = 1;
for i=1:length(dset) % for each pedestrian
	
	BlobsN(ii).mvec	= [];
	BlobsN(ii).pvec	= [];
	BlobsN(ii).index= [];
	BlobsN(ii).mvec	= Blobs(dset(i).globalindex(1)).mvec; %
	BlobsN(ii).pvec	= Blobs(dset(i).globalindex(1)).pvec;
	BlobsN(ii).index = dset(i).globalindex(1);
	ii = ii+1;
	BlobsN(ii).mvec	= Blobs(dset(i).globalindex(end)).mvec; % static feature
	BlobsN(ii).pvec	= Blobs(dset(i).globalindex(end)).pvec;
	BlobsN(ii).index	= dset(i).globalindex(end);
	ii = ii+1;
end


%% MSCR distances computation
hh = waitbar(0,'Image features computation...');
for i=1:length(BlobsN)
    %%---- MSE analysis
    mser    =   BlobsN(i);

    [keep_list,BB]= blobs_inside_image(mser.mvec,[H,W]); % used to compute BB
    keep_list = [1:size(mser.mvec,2)];
	
	centroid = zeros(2,length(keep_list));
	colour = zeros(3,length(keep_list));
    for j=1:length(keep_list);
        centroid(:,j) = BB(1:2,1,j);
        colour(:,j) = applycform(mser.pvec(:,j)', reg);
	end

	BlobsN(i).Mmvec=mser.mvec(:,keep_list);
    BlobsN(i).Mpvec=colour;
	
	% blob clustering 
% 	data = [centroid;colour];
% 	[bestk,bestpp,bestmu,bestcov,dl,countf] = mixtures4(data,kmin,min(kmax,size([centroid;colour],2)),regularize,th,covoption);
% 	data = data';
% 	num_data = size(data,1);
% 	label	 = zeros(num_data,1);
% 	ll		 = zeros(num_data,1);
% 	for j = 1:num_data
% 		vote = zeros(bestk,1);
% 		for k=1:bestk
% 			vote(k)= mvnpdf(data(j,:)',bestmu(:,k),bestcov(:,:,k));
% 		end
% 		[ll(j), label(j)] = max(vote);
% 	end
% 	% 1-NN
% 	imax = zeros(bestk,1);
% 	for k=1:bestk
% 		tmp			= zeros(length(ll),1);
% 		inds		= (label == k);
% 		tmp(inds)	= ll(inds);
% 		[trashm, imax(k)] = max(tmp);
% 	end
%     BlobsN(i).Mmvec=mser.mvec(:,imax);
%     BlobsN(i).Mpvec=colour(:,imax);
		
% 	
% 	mser    =   BlobsS(i);
% 	
% 	[keep_list,BB]= blobs_inside_image(mser.mvec,[H,W]); % used to compute BB
% 	keep_list = [1:size(mser.mvec,2)];
% 	
% 	centroid = zeros(2,length(keep_list));
% 	colour = zeros(3,length(keep_list));
% 	for j=1:length(keep_list);
% 		centroid(:,j) = BB(1:2,1,j);
% 		colour(:,j) = applycform(mser.pvec(:,j)', reg);
%     end
% 	
% 	BlobsS(i).Mmvec=mser.mvec(:,keep_list);
%     BlobsS(i).Mpvec=colour;

    waitbar(i/length(dset),hh);
	
end
close(hh);

%%
hh = waitbar(0,'Distances computation MSCR...');
for i = 1:length(BlobsN)
	
	Mmvec{1} = BlobsN(i).Mmvec(3,:);
	Mpvec{1} = BlobsN(i).Mpvec;
	
	num(1)=size(Mmvec{1},2);
	
	dist_y_n	= cell(length(dset),1);
	dist_color_n= cell(length(dset),1);
	
	for j = 1:length(BlobsN)
		
		if i~=j
			Mmvec{2} = BlobsN(j).Mmvec(3,:);
			Mpvec{2} = BlobsN(j).Mpvec;
			
			num(2)= size(Mmvec{2},2); %length(Mindeces{2});
			[thrash,max_ind(i,j)]= max([num(1),num(2)]);
			max_i = mod(max_ind(i,j),2)+1;%max_ind;
			min_i = max_ind(i,j);%mod(max_ind,2)+1;
			max_info = num(max_i);
			min_info = num(min_i);
			
			dist_y{j} = 1000*ones(min_info, max_info);
			dist_color{j} = 1000*ones(min_info, max_info);
			
			% distanze tutti i blob di i con tutti i blob di j
			for k=1:max_info
				for h=1:min_info
					
					dist_y{j}(h,k) = abs(Mmvec{max_i}(k)-Mmvec{min_i}(h));
					dist_color{j}(h,k) = sqrt(sum((Mpvec{max_i}(:,k)-Mpvec{min_i}(:,h)).^2));
					
				end
			end
			
			
			%% check outliers
			
			ref_y = min(dist_y{j}); me_ref_y = mean(ref_y); std_ref_y = std(ref_y);
			ref_color = min(dist_color{j}); me_ref_color = mean(ref_color); std_ref_color = std(ref_color);
			good = find((ref_y<(me_ref_y+3.5*std_ref_y))&(ref_color<(me_ref_color+3.5*std_ref_color)));
			
			if ~isempty(good)
				max_useful_info = length(good);
				dist_y2 = dist_y{j}(:,good);
				dist_color2 = dist_color{j}(:,good);
				
				%%normalize
				dist_y_n{j} = dist_y{j}./max(dist_y2(:));
				dist_color_n{j} = dist_color{j}./max(dist_color2(:));
				
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
		else
			final_dist_y(i,j)		= 0;
			final_dist_color(i,j)	= 0;
		end
	end
	waitbar(i*j/length(dset)^2,hh);
	clear dist_y;
	clear dist_color;
	clear dist_color_n;
	clear dist_y_n;
	
end
close(hh);


for i=1:length(dset)
    final_dist_y(i,:) =  final_dist_y(i,:)./max(final_dist_y(i,:));
    final_dist_color(i,:) = final_dist_color(i,:)./max(final_dist_color(i,:));
    final_mscr_dist(i,:) = (pyy*final_dist_y(i,:) +  pcc*final_dist_color(i,:) ); % +  final_dist_hist(i,:); 
end

%% 10 times validation
for t = 1:10

	S1 = 1:length(dset); % all data (1 vs. 1, dyn vs. static)
	S2 = length(dset)+1:length(dset)+length(S1);
	
    final_dist_y_tmp    = final_dist_y(S1,S2);
    final_dist_color_tmp= final_dist_color(S1,S2);
%     final_dist_hist_tmp = final_dist_hist(S1,S2);   
%     final_dist_epitext  = dist_epitext(S1,S2);
    
    for i=1:length(S1)
        final_dist_y_tmp(i,:) =  final_dist_y_tmp(i,:)./max(final_dist_y_tmp(i,:));
        final_dist_color_tmp(i,:) =  final_dist_color_tmp(i,:)./max(final_dist_color_tmp(i,:));
        final_dist_tmp(i,:) = (pyy*final_dist_y_tmp(i,:) + pcc*final_dist_color_tmp(i,:)); % +  final_dist_hist_tmp(i,:) + pee*final_dist_epitext(i,:);
        [dists, ordered] = sort(final_dist_tmp(i,:),'ascend');
        punt_final_dist_tmp(i)= find(ordered==i);
    end

    %% CMC
    for i=1:length(S1)
        CMC(i) = length(find(punt_final_dist_tmp<=i));
    end
    AUC = sum(CMC);
    nAUC = sum(CMC)/(length(S1)*length(S1))*100;

    %% SRR
    M = 1:10;
    for m = 1:length(M)
        SRR(1,m) = CMC(uint16(floor(length(S1)/M(m))))/length(S1)*100;
    end

    stats.A1(t,:)    = [S1];
    stats.A2(t,:)    = [S2];
    stats.CMC(t,:)      = CMC;
    stats.AUC(t)        = AUC;
    stats.nAUC(t)       = nAUC;
    stats.SRR(t,:)      = SRR;
    stats.ordermatch(t,:)    = punt_final_dist_tmp;
end

fprintf('AUC: %d     normalized AUC: %f \n',mean(stats.AUC), mean(stats.nAUC))
figure, plot(mean(stats.CMC,1)/length(S1)*100, 'b.-', 'Linewidth',1.2); title('Cumulative Matching Characteristic (CMC)'), grid on, axis([1,length(S1),0,100]);
figure,plot(M,mean(stats.SRR,1),'o-'),title('Synthetic Recognition Rate (SRR)'), grid on, axis([min(M),max(M),0,100])
