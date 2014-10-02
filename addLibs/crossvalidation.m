%% validation for dynamic feature
clear stats 

switch dataname
	
	case 'iLIDS' %% select 2 view 10 times, randomply
		for t = 1:10
			
			[S1,S2] = ExtractDataset(dset);
            
            for i=1:length(S1)
                n = dir_e_list(permit_inds(S1(i))).name;
                id = str2num(n(1:4));
                idindex(i) = id;
            end
            
            
            %S2 = setdiff([1:length(permit_inds)],S1');
            n = dir_e_list(permit_inds(S2(end))).name;
            idpersonmax = str2num(n(1:4));
            
            for i=1:idpersonmax
                Score(i).match = [];
            end
            
            for k =1:length(S2)
                
                n = dir_e_list(permit_inds(S2(k))).name;
                idperson = str2num(n(1:4));
                pos = find(idindex == idperson);
                
                final_dist_y_tmp    = final_dist_y(S2(k),S1);
                final_dist_color_tmp= final_dist_color(S2(k),S1);
                final_dist_hist_tmp = final_dist_hist(S2(k),S1);
                final_dist_epitext_tmp  = dist_epitext(S2(k),S1);
                final_dist_tmp = (pyy*final_dist_y_tmp + pcc*final_dist_color_tmp) +  final_dist_hist_tmp + pee*final_dist_epitext_tmp;
                [dists, ordered] = sort(final_dist_tmp,'ascend');
                
                a = find(ordered==pos);
                if  isempty(Score(idperson).match)
                    Score(idperson).match = a;
                else
                    Score(idperson).match = [Score(idperson).match a];
                end
            end
            
           
            j=1;
            for k=1:idpersonmax
                if ~isempty(Score(k).match)
                    for m=1:length(Score(k).match)
                        punt_final_dist_tmp(j) = Score(k).match(m);
                        j = j+1;
                    end
                end
            end
            
			
			%% CMC
			CMC = zeros(length(S1),1);
			for i=1:length(S1)
				CMC(i) = length(find(punt_final_dist_tmp<=i));
			end
			AUC = sum(CMC);
			nAUC = sum(CMC)/(length(S1)*length(S2))*100;
			
			%% SRR
			M = 1:10;
			SRR = zeros(1,length(M));
			for m = 1:length(M)
				SRR(1,m) = CMC(uint16(floor(length(S1)/M(m))))/length(S2)*100;
			end
			
			stats.A1(t,:)    = [S1];
			%stats.A2(t,:)    = [S2];
			stats.CMC(t,:)      = CMC;
			stats.AUC(t)        = AUC;
			stats.nAUC(t)       = nAUC;
			stats.SRR(t,:)      = SRR;
		end
		
				
	case 'VIPeR' %% 5x2 cross-validation: split in two the dataset and use only one partition
		for t = 1:10

            clear final_dist_tmp punt_final_dist_tmp CMC
            
            [S1,S2] = ExtractDatasetOrd(dset);

            indimg = randperm(length(cvpridx))';
            indimg = cvpridx(indimg); 
            indimg = indimg(1:316); % compare against Zeng et al. (group re-id)
%             indimg = indimg(1:474); % compare against Prosser et al. (support vector ranking)
            
            S1 = S1(indimg);
            S2 = S2(indimg);
			
			final_dist_y_tmp    = final_dist_y(S1,S2);
			final_dist_color_tmp= final_dist_color(S1,S2);
			final_dist_hist_tmp = final_dist_hist(S1,S2);
 			final_dist_epitext  = dist_epitext(S1,S2);
			
            
            for i=1:length(S1)
				final_dist_y_tmp(i,:) =  final_dist_y_tmp(i,:)./max(final_dist_y_tmp(i,:));
				final_dist_color_tmp(i,:) =  final_dist_color_tmp(i,:)./max(final_dist_color_tmp(i,:));
				final_dist_tmp(i,:) = (pyy*final_dist_y_tmp(i,:) + pcc*final_dist_color_tmp(i,:)) +  final_dist_hist_tmp(i,:)+ pee*final_dist_epitext(i,:);
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
			M = 1:25;
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
			
	case {'ETHZ1','ETHZ2','ETHZ3'} %% select the best view for first cam, and the others randomly
		for t = 1:100
			
            [S1,bleah] = ExtractDataset(dset);
            
            for i=1:length(S1)
                n = dir_e_list(permit_inds(S1(i))).name;
                id = str2num(n(1:4));
                idindex(i) = id;
            end
            
            
            S2 = setdiff([1:length(permit_inds)],S1');
            n = dir_e_list(permit_inds(S2(end))).name;
            idpersonmax = str2num(n(1:4));
            
            for i=1:idpersonmax
                Score(i).match = [];
            end
            
            for k =1:length(S2)
                
                n = dir_e_list(permit_inds(S2(k))).name;
                idperson = str2num(n(1:4));
                pos = find(idindex == idperson);
                
                final_dist_y_tmp    = final_dist_y(S2(k),S1);
                final_dist_color_tmp= final_dist_color(S2(k),S1);
                final_dist_hist_tmp = final_dist_hist(S2(k),S1);
                final_dist_epitext_tmp  = dist_epitext(S2(k),S1);
                final_dist_tmp = (pyy*final_dist_y_tmp + pcc*final_dist_color_tmp) +  final_dist_hist_tmp + pee*final_dist_epitext_tmp;
                final_dist_visu(k,:) = final_dist_tmp;
                [dists, ordered] = sort(final_dist_tmp,'ascend');
                
                a = find(ordered==pos);
                if  isempty(Score(idperson).match)
                    Score(idperson).match = a;
                else
                    Score(idperson).match = [Score(idperson).match a];
                end
            end
            
            %     j = 1;
            %     for k=1:idpersonmax
            %         if ~isempty(Score(k).match)
            %             punt_final_dist_tmp(j) = min(Score(k).match);
            %             j = j+1;
            %         end
            %     end
            
            j=1;
            for k=1:idpersonmax
                if ~isempty(Score(k).match)
                    for m=1:length(Score(k).match)
                        punt_final_dist_tmp(j) = Score(k).match(m);
                        j = j+1;
                    end
                end
            end
            

			
			%% CMC
			CMC = zeros(length(S1),1);
			for i=1:length(S1)
				CMC(i) = length(find(punt_final_dist_tmp<=i));
			end
			AUC = sum(CMC);
			nAUC = sum(CMC)/(length(S1)*length(S2))*100;
			
			%% SRR
			M = 1:10;
			SRR = zeros(1,length(M));
			for m = 1:length(M)
				SRR(1,m) = CMC(uint16(floor(length(S1)/M(m))))/length(S2)*100;
			end
			
			stats.A1(t,:)    = [S1];
			%stats.A2(t,:)    = [S2];
			stats.CMC(t,:)      = CMC;
			stats.AUC(t)        = AUC;
			stats.nAUC(t)       = nAUC;
			stats.SRR(t,:)      = SRR;
			%stats.ordermatch(t,:)    = punt_final_dist_tmp;
		end
end

fprintf('AUC: %d     normalized AUC: %f \n',mean(stats.AUC), mean(stats.nAUC))
figure(1), plot(mean(stats.CMC,1)/length(S2)*100, 'b.-', 'Linewidth',1.2); title('Cumulative Matching Characteristic (CMC)'), grid on, axis([1,50,0,100]);
%figure(2),plot(M,mean(stats.SRR,1),'o-'),title('Synthetic Recognition Rate (SRR)'), grid on, axis([min(M),max(M),0,100])

            
            
% 			S1 = find(selected_models == permits_ind); %% CHECK THIS!!!!
% 			
% 			%S2 = 2:2:length(dset)*2; % even
% 			[Bleah,S2] = ExtractDataset(dset);
% 			
% 			
% 			final_dist_tmp = zeros(length(S1),length(S2));
% 			punt_final_dist_tmp = zeros(length(S1),1);
% 			
% 			final_dist_y_tmp    = final_dist_y(S1,S2);
% 			final_dist_color_tmp= final_dist_color(S1,S2);
% 			final_dist_hist_tmp = final_dist_hist(S1,S2);
% 			%     final_dist_epitext  = dist_epitext(S1,S2);
% 			
% 			for i=1:length(S1)
% 				final_dist_y_tmp(i,:) =  final_dist_y_tmp(i,:)./max(final_dist_y_tmp(i,:));
% 				final_dist_color_tmp(i,:) =  final_dist_color_tmp(i,:)./max(final_dist_color_tmp(i,:));
% 				final_dist_tmp(i,:) = (pyy*final_dist_y_tmp(i,:) + pcc*final_dist_color_tmp(i,:)) +  final_dist_hist_tmp(i,:);% + pee*final_dist_epitext(i,:);
% 				[dists, ordered] = sort(final_dist_tmp(i,:),'ascend');
% 				punt_final_dist_tmp(i)= find(ordered==i);
% 			end
% 			
% 			%% CMC
% 			CMC = zeros(length(S1),1);
% 			for i=1:length(S1)
% 				CMC(i) = length(find(punt_final_dist_tmp<=i));
% 			end
% 			AUC = sum(CMC);
% 			nAUC = sum(CMC)/(length(S1)*length(S1))*100;
% 			
% 			%% SRR
% 			M = 1:10;
% 			SRR = zeros(1,length(M));
% 			for m = 1:length(M)
% 				SRR(1,m) = CMC(uint16(floor(length(S1)/M(m))))/length(S1)*100;
% 			end
% 			
% 			stats.A1(t,:)    = [S1];
% 			stats.A2(t,:)    = [S2];
% 			stats.CMC(t,:)      = CMC;
% 			stats.AUC(t)        = AUC;
% 			stats.nAUC(t)       = nAUC;
% 			stats.SRR(t,:)      = SRR;
% 			stats.ordermatch(t,:)    = punt_final_dist_tmp;
% 		end
% end
% 
% fprintf('AUC: %d     normalized AUC: %f \n',mean(stats.AUC), mean(stats.nAUC))
% figure, plot(mean(stats.CMC,1)/size(stats.CMC,2)*100, 'b.-', 'Linewidth',1.2); title('Cumulative Matching Characteristic (CMC)'), grid on, axis([1,size(stats.CMC,2),0,100]);
% figure,plot(M,mean(stats.SRR,1),'o-'),title('Synthetic Recognition Rate (SRR)'), grid on, axis([min(M),max(M),0,100])
