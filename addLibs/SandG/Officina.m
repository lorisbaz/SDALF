%% Officina 20_10_2009
%% ----- Con learning gaussiano caso modello statico -- modello ------
clear all, clc
load TMP_gauss_clust

avoid_inds  = [11,14,33,42,60,67,68,91,92,96,127,128,147,148,149,165,195,210,...
	212,217,218,236,248,259,260,290,303,310,328,374,427,439,39,110,120,166,...
	167,178,338];

dir_e   =   '/Users/bazzani/Pictures/i-LIDS_Pedestrian'; % '/Volumes/Untitled 1/datasets/Tagging/i-LIDS_Pedestrian';

%% Init
dir_e_list = dir(strcat(dir_e,'/*.jpg'));
n_img       =   size(dir_e_list,1);
permit_inds = setdiff(1:n_img,avoid_inds);

kmin=1;
kmax=5;
regularize = 0;
th = 1e-4;
covoption=1;


%% TRAINING: automatic extraction of outlier threashold
for p = 1:3
	for i = 1:length(permit_inds) % BY NOW, we use all the dataset (not in the future!)
		y = [MSCR(i,p).data];
		if ~isempty(y)
			[bestk,bestpp,bestmu,bestcov,dl,countf] = mixtures4(y,kmin,min(kmax,size(y,2)),regularize,th,covoption);
			model{i,p}.cov		= bestcov;
			model{i,p}.mu		= bestmu;
			model{i,p}.alpha	= bestpp;
			model{i,p}.K		= bestk;
			% controlla la dimensionalit? degli output
			data = y';
			num_data = size(data,1);
			label	 = zeros(num_data,1);
			for j = 1:num_data
				vote = zeros(bestk,1);
				for k=1:bestk
					vote(k)= mvnpdf(data(j,:)',bestmu(:,k),bestcov(:,:,k));
				end
				[trash, label(j)] = max(vote);
			end
			model{i,p}.label	= label;
			cov(:,:,i) = mean(bestcov,3); % automatic extraction of outlier threashold
		else
			model{i,p} = [];
			cov(:,:,i) = zeros(4,4);
		end
	end
	meancov(:,:,p) = mean(cov,3);
end

%% TESTING
hh = waitbar(0,'Distances computation MSCR...');
for p = 1:3
	for i = 1:length(permit_inds)
		if ~isempty(model{i,p})
			data = [MSCR(i,p).data]';
			for j = 1:length(permit_inds)
				tdata = [MSCR(j,p).data]';
				if ~isempty(tdata)
				num_tdata = size(tdata,1);
				label	 = zeros(num_tdata,1);
				% associate each blob to one cluster
				for n = 1:num_tdata
					vote = zeros(model{i}.K,1);
					for k=1:model{i}.K
						% outliers elimination using learned covariance matrix
						vote2= mvnpdf(tdata(n,:)',model{i,p}.mu(:,k),meancov(:,:,p)); 
						if (vote2<mvnpdf(diag(2*sqrt(meancov(:,:,p))),model{i,p}.mu(:,k),meancov(:,:,p)))
							vote(k) = 0;
						else % no outlier: voting
							vote(k) = mvnpdf(tdata(n,:)',model{i,p}.mu(:,k),model{i,p}.cov(:,:,k));
						end
					end
					[trash, label(n)] = max(vote);
					
					if label(n) ~=0
						tmp = data(model{i,p}.label == label(n),:);
						for k=1:size(tmp,1)
							dist_y{j}(n,k) = abs(tmp(k,1)-tdata(n,2));
							dist_color{j}(n,k) = sqrt(sum((tmp(k,2:4)-tdata(n,2:4)).^2));
						end
					end
				end
				
				final_dist_yp(i,j,p)  = sum(dist_y{j})/max(dist_y{j});
				final_dist_colorp(i,j,p) = sum(dist_color{j})/max(dist_color{j});
				
				waitbar((i*j*p)/length(permit_inds)^2*3,hh);
				else
					final_dist_yp(i,j,p)  = 1;
					final_dist_colorp(i,j,p)  = 1;
				end
			end
		else
			final_dist_yp(i,:,p) = ones(1,length(permit_inds));
			final_dist_colorp(i,:,p) = ones(1,length(permit_inds));
		end
	end
	final_dist_y(i,j) = sum(final_dist_yp(i,j,p),3)/3;
	final_dist_color(i,j) = sum(final_dist_colorp(i,j,p),3)/3;
end
close(hh);

%% 10 times validation
SetDataset;
for t = 1:10

    [S1,S2] = ExtractDataset(dset);
    
    final_dist_y_tmp    = final_dist_y(S1,S2);
    final_dist_color_tmp= final_dist_color(S1,S2);
    
    for i=1:length(S1)
        final_dist_y_tmp(i,:) =  final_dist_y_tmp(i,:)./max(final_dist_y_tmp(i,:));
        final_dist_color_tmp(i,:) =  final_dist_color_tmp(i,:)./max(final_dist_color_tmp(i,:));
        final_dist_tmp(i,:) = (pyy*final_dist_y_tmp(i,:) + pcc*final_dist_color_tmp(i,:));
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

