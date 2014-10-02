%% validation for dynamic feature (dynamic feature computing is done before)


final_dist_tmp = zeros(length(statmodels),length(dset));
punt_final_dist_tmp = zeros(length(statmodels),1);

for i=1:length(statmodels)
    
    n = dir_e_list(statmodels(i)).name;
    idperson = str2num(n(1:4));
   
    
	final_dist_y_tmp(i,:) =  final_dist_y(i,:)./max(final_dist_y(i,:));
	final_dist_color_tmp(i,:) =  final_dist_color(i,:)./max(final_dist_color(i,:));
	final_dist_tmp(i,:) = (pyy*final_dist_y_tmp(i,:) + pcc*final_dist_color_tmp(i,:)) +  final_dist_hist(i,:) + pee * dist_epitext(i,:);  % if you want epiteture decomment this
	[dists, ordered] = sort(final_dist_tmp(i,:),'ascend');
    
	punt_final_dist_tmp(i)= find(ordered==idperson);
end

%% CMC
CMC = zeros(length(dset),1);
for i=1:length(dset)
	CMC(i) = length(find(punt_final_dist_tmp<=i));
end
AUC = sum(CMC);
nAUC = sum(CMC)/(length(dset)*length(statmodels))*100;

%% SRR
M = 1:min(10,length(dset));
SRR = zeros(1,length(M));
for m = 1:length(M)
	SRR(1,m) = CMC(uint16(floor(length(dset)/M(m))))/length(statmodels)*100;
end

stats.A1{t}    = [dynmodels];
stats.CMCn(t,:) = CMC/length(statmodels) * 100;
stats.CMC(t,:)      = CMC;
stats.AUC(t)        = AUC;
stats.nAUC(t)       = nAUC;
stats.SRR(t,:)      = SRR;
stats.ordermatch(t,:)    = punt_final_dist_tmp;

