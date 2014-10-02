%% maximally-texturized distances computation 
%load(['TxPatch_f' num2str(SUBfac) '.mat']);
hwait	= waitbar(0,'Distance between epitextures...');
part	= 1; % upper-body part
for i=1:length(permit_inds)
   
  if ~isempty(max_txpatch(i,part).lbph) && ~isempty(max_txpatch(i,part).lbph{1} ) 
    
      for j=1:length(permit_inds)
        
        db = []; dl = [];
        
        if ~isempty( max_txpatch(j,part).lbph)  && ~isempty( max_txpatch(j,part).lbph{1})  
            
            for h=1:length(max_txpatch(i,part).lbph)
                for k=1:length(max_txpatch(j,part).lbph)    
                    db(h,k)= bhattacharyya(max_txpatch(i,part).lbph{h}, max_txpatch(j,part).lbph{k});
                end
            end
        end    
        
        
        if ~isempty(db) 
            dist_txpatch(i,j) = min(db(:));
        else
            dist_txpatch(i,j) = 0.5;
        end
	  end
  else
	  dist_txpatch(i,1:length(permit_inds)) = 0.5; 
  end
  
  waitbar(i/length(permit_inds),hwait);
end

close(hwait);

dist_epitext = dist_txpatch;  

% name = ['iLIDS_TxPatchmatch_f' num2str(SUBfac) '.mat'];
% save(name, 'dist_txpatch');
