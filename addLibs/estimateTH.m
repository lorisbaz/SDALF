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

a = diag(final_dist_tmp);
b = [final_dist_tmp-diag(a) ];
b = vec(b(find(b~=0)));


x = (0:0.01:3)';
y1 = gaussmf(x, [std(a) mean(a)]);
y2 = gaussmf(x, [std(b) mean(b)]);

plot(x,[y1,y2]);

