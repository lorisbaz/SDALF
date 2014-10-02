hh = waitbar(0,'Distances computation RGB images...');

for i = 1:n_img
    Mmvec{1} = dataf(i,1).Mmvec(3,:);
    Mpvec{1} = dataf(i,1).Mpvec;
    num(1)=size(Mmvec{1},2);

    for j = 1:n_img
        
        Mmvec{2} = dataf(j,2).Mmvec(3,:);
        Mpvec{2} = dataf(j,2).Mpvec;
        num(2)= size(Mmvec{2},2); %length(Mindeces{2});
        
        [thrash,max_ind(i,j)]= max([num(1),num(2)]);
        max_i = mod(max_ind(i,j),2)+1;%max_ind;
        min_i = max_ind(i,j);%mod(max_ind,2)+1;
        max_info = num(max_i);
        min_info = num(min_i);
       
        dist_y = []; 
        dist_color = [];
       
        
        % distanze tutti con tutti
        for k=1:max_info           
              for h=1:min_info
                  dist_y(h,k) = abs(Mmvec{max_i}(k)-Mmvec{min_i}(h));
                  dist_color(h,k) = sqrt(sum((Mpvec{max_i}(:,k)-Mpvec{min_i}(:,h)).^2));
                  
              end
        end
        
        
        %% check outliers 
        
        ref_y = min(dist_y); me_ref_y = mean(ref_y); std_ref_y = std(ref_y);
        ref_color = min(dist_color); me_ref_color = mean(ref_color); std_ref_color = std(ref_color);
        good = find((ref_y<(me_ref_y+2.5*std_ref_y)) &(ref_color<(me_ref_color+2.5*std_ref_color)));
        max_useful_info = length(good);

        dist_y2 = dist_y(:,good);
        dist_color2 = dist_color(:,good);
        
        %%normalize
        dist_y_n = (dist_y)./max(dist_y2(:));
        dist_color_n = (dist_color)./max(dist_color2(:));
     
        %%distance computation
        totdist_n = ( pyy*dist_y_n(:,good)  + pcc*dist_color_n(:,good) ) ;
        
        %%Minimization
        [min_dist,matching] = min(totdist_n);
        
        useful_i = sub2ind(size(totdist_n),[matching]',[1:max_useful_info]'); 
        final_dist_y(i,j)  = sum(dist_y2(useful_i))/max_useful_info;
        final_dist_color(i,j) = sum(dist_color2(useful_i))/max_useful_info;
        match(i,j).info = [matching; good];
        match(i,j).dist = min_dist;
        
        %final_dist_hist(i,j) = bhattacharyya(dataf(i,1).Hist, dataf(j,2).Hist);
    end
    waitbar(i/n_img,hh);
    clear dist_y;     clear dist_color;
    clear dist_color_n;
    clear dist_y_n; 

end