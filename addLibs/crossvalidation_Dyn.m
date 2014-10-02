%% validation for dynamic feature (dynamic feature computing is done before)
S1 = 1:2:length(dset)*2; % odd
S2 = 2:2:length(dset)*2; % even

final_dist_tmp = zeros(length(S1),length(S2));
punt_final_dist_tmp = zeros(length(S1),1);

final_dist_y_tmp    = final_dist_y(S1,S2);
final_dist_color_tmp= final_dist_color(S1,S2);
final_dist_hist_tmp = final_dist_hist(S1,S2);
final_dist_epitext  = dist_epitext(S1,S2);  % if you want epiteture decomment this

for i=1:length(S1)
	final_dist_y_tmp(i,:) =  final_dist_y_tmp(i,:)./max(final_dist_y_tmp(i,:));
	final_dist_color_tmp(i,:) =  final_dist_color_tmp(i,:)./max(final_dist_color_tmp(i,:));
	final_dist_tmp(i,:) = (pyy*final_dist_y_tmp(i,:) + pcc*final_dist_color_tmp(i,:)) +  final_dist_hist_tmp(i,:) + pee * final_dist_epitext(i,:);  % if you want epiteture decomment this
	[dists, ordered] = sort(final_dist_tmp(i,:),'ascend');
	punt_final_dist_tmp(i)= find(ordered==i);
end

%% CMC
CMC = zeros(length(S1),1);
for i=1:length(S1)
	CMC(i) = length(find(punt_final_dist_tmp<=i));
end
AUC = sum(CMC);
nAUC = sum(CMC)/(length(S1)*length(S1))*100;

%% SRR
M = 1:min(10,length(S1));
SRR = zeros(1,length(M));
for m = 1:length(M)
	SRR(1,m) = CMC(uint16(floor(length(S1)/M(m))))/length(S1)*100;
end

stats.A1(t,:)    = [S1];
stats.A2(t,:)    = [S2];
stats.dists(t,:,:)  = final_dist_tmp;
stats.CMC(t,:)      = CMC;
stats.AUC(t)        = AUC;
stats.nAUC(t)       = nAUC;
stats.SRR(t,:)      = SRR;
stats.ordermatch(t,:)    = punt_final_dist_tmp;

