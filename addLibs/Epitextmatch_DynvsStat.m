%% Faccio clustering delle patch appartenenti alle immagini delle signature
clear max_txpatchDyn
clear dist_txpatchDyn

ii = 1; ddtmp = 1;
for i = 1:length(dset) % for each ped
    
    % -- dyn 1
    max_txpatchDyn(ii).lbph = [];
    dd = length(ped(i).dynmodels) - 1;
    rnddd = ped(i).rnddd;
    ep = 1;
    for j = 1:dd % for each view/2 
        ind = find(ped(i).dynmodels(rnddd(j)) == dset(i).globalindex); % we are searching for the local index (in permits_ind)
        if ~isempty(max_txpatch(dset(i).index(ind),2))
            for m = 1:length(max_txpatch(dset(i).index(ind),2).lbph)
                max_txpatchDyn(ii).lbph{ep}  = max_txpatch(dset(i).index(ind),2).lbph{m};
                ep = ep+1;
            end
        end
    end
    ii = ii+1;
    
     % -- dyn 2
    max_txpatchDyn(ii).lbph = [];
    %dd = uint16(length(ped(i).dynmodels)/2);
    %rnddd = ped(i).rnddd;
    ep = 1;
    for j = dd+1:length(ped(i).dynmodels) % for each view/2
        ind = find(ped(i).dynmodels(rnddd(j)) == dset(i).globalindex); % we are searching for the local index (in permits_ind)
        if ~isempty(max_txpatch(dset(i).index(ind),2))
            for m = 1:length(max_txpatch(dset(i).index(ind),2).lbph)
                max_txpatchDyn(ii).lbph{ep}  = max_txpatch(dset(i).index(ind),2).lbph{m};
                ep = ep+1;
            end
        end
    end
    ii = ii+1;
    
    
    
end


% maximally-texturized distances computation 
%load(['TxPatch_f' num2str(SUBfac) '.mat']);
hwait	= waitbar(0,'Distance between epitextures...');
part	= 1; % upper-body part
for i=1:size(max_txpatchDyn,2)
   
  if ~isempty(max_txpatchDyn(i).lbph) && ~isempty(max_txpatchDyn(i).lbph{1} ) 
    
      for j=1:size(max_txpatchDyn,2)
        
        db = []; dl = [];
        
        if ~isempty( max_txpatchDyn(j).lbph)  && ~isempty( max_txpatchDyn(j).lbph{1})  
            
            for h=1:length(max_txpatchDyn(i).lbph)
                for k=1:length(max_txpatchDyn(j).lbph)    
                    db(h,k)= bhattacharyya(max_txpatchDyn(i).lbph{h}, max_txpatchDyn(j).lbph{k});
                end
            end
        end    
        
        
        if ~isempty(db) 
            dist_txpatchDyn(i,j) = min(db(:));
        else
            dist_txpatchDyn(i,j) = 0.5;
        end
	  end
  else
	  dist_txpatchDyn(i,1:size(max_txpatchDyn,2)) = 0.5; 
  end
  
  waitbar(i/size(max_txpatchDyn,2),hwait);
end

close(hwait);

dist_epitext = dist_txpatchDyn;  

% name = ['iLIDS_TxPatchmatch_f' num2str(SUBfac) '.mat'];
% save(name, 'dist_txpatch');



% %% Faccio clustering delle patch appartenenti alle immagini delle signature
% 
% ii = 1; ddtmp = 1;
% for i = 1:length(dset) % for each ped
%     
%     max_txpatchDyn(ii).lbph = [];
%     dd = uint16(length(ped(i).dynmodels));
%     rnddd = ped(i).rnddd;
%     ep = 1;
%     for j = 1:dd % for each view/2 
%         ind = find(ped(i).dynmodels(rnddd(j)) == dset(i).globalindex); % we are searching for the local index (in permits_ind)
%         if ~isempty(max_txpatch(dset(i).index(ind),2))
%             for m = 1:length(max_txpatch(dset(i).index(ind),2).lbph)
%                 max_txpatchDyn(ii).lbph{ep}  = max_txpatch(dset(i).index(ind),2).lbph{m};
%                 ep = ep+1;
%             end
%         end
%     end
%     ii = ii+1;
% end   
% 
% for i=1:length(statmodels)
%     max_txpatchStat(i).lbph = [];
%     ep = 1;
%     if ~isempty(max_txpatch(statmodels(i),2))
%         for m = 1:length(max_txpatch(statmodels(i),2).lbph)
%             max_txpatchStat(i).lbph{ep}  = max_txpatch(statmodels(i),2).lbph{m};
%             ep = ep+1;
%         end
%     end
% end
% 
% % maximally-texturized distances computation 
% %load(['TxPatch_f' num2str(SUBfac) '.mat']);
% hwait	= waitbar(0,'Distance between epitextures...');
% part	= 1; % upper-body part
% for i=1:size(max_txpatchStat,2)
%    
%   if ~isempty(max_txpatchStat(i).lbph) && ~isempty(max_txpatchStat(i).lbph{1} ) 
%     
%       for j=1:size(max_txpatchDyn,2)
%         
%         db = []; dl = [];
%         
%         if ~isempty( max_txpatchDyn(j).lbph)  && ~isempty( max_txpatchDyn(j).lbph{1})  
%             
%             for h=1:length(max_txpatchStat(i).lbph)
%                 for k=1:length(max_txpatchDyn(j).lbph)    
%                     db(h,k)= bhattacharyya(max_txpatchStat(i).lbph{h}, max_txpatchDyn(j).lbph{k});
%                 end
%             end
%         end    
%         
%         
%         if ~isempty(db) 
%             dist_txpatchDyn(i,j) = min(db(:));
%         else
%             dist_txpatchDyn(i,j) = 0.5;
%         end
% 	  end
%   else
% 	  dist_txpatchDyn(i,1:size(max_txpatchDyn,2)) = 0.5; 
%   end
%   
%   waitbar(i/size(max_txpatchStat,2),hwait);
% end
% 
% close(hwait);
% 
% dist_epitext = dist_txpatchDyn;  
% 
% % name = ['iLIDS_TxPatchmatch_f' num2str(SUBfac) '.mat'];
% % save(name, 'dist_txpatch');
