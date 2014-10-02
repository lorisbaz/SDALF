%% Extract epitexture
Nep = 20;
Nth = 4;
patchSize = [21*SUBfac,11*SUBfac]; % default
% mapping     = getmapping(8,'u2'); % lbp
switch dataname
    
    case 'ETHZ1'
        load(['MAT/Epitome_ETZH1' '_f' num2str(SUBfac) '_Exp' expnum '.mat']);
        hh = waitbar(0,'Extract Epitoms...');
		for i = 1:size(mu,4)
			ordPatches = likelyPatches( mu(:,:,:,i), pcl(:,:,i), patchSize, Nep );
			histo_entropy = zeros(size(ordPatches,4),1);
			for p = 1:size(ordPatches,4)
				histo_entropy(p)  = entropy(ordPatches(:,:,1,p)) + entropy(ordPatches(:,:,2,p)) + entropy(ordPatches(:,:,3,p));
				% 		histo_entropy(p)  = entropy(rgb2gray(ordPatches(:,:,:,p)));
			end
			[dum ord] = sort(histo_entropy,'descend');
			ord = ord(1:Nth);
			epitexture(i).patch(:,:,:,1:Nth) = ordPatches(:,:,:,ord); % even
			for p = 1:Nth
				img_hsv     =   rgb2hsv(ordPatches(:,:,:,ord(p)));
				epitexture(i).hist(:,p)=   [hist(raster(img_hsv(:,:,1)),[0:1/NBINs(1):1]) ...
					hist(raster(img_hsv(:,:,2)),[0:1/NBINs(2):1])...
					hist(raster(img_hsv(:,:,3)),[0:1/NBINs(3):1])];
% 				epitexture(i).lbph(:,p) = lbp(rgb2gray(ordPatches(:,:,:,ord(p))), 1,8,mapping, 'h');
			end
			waitbar(i/(size(mu,4)),hh);
		end
		close(hh);
        
    case 'ETHZ2'
        load(['MAT/Epitome_ETHZ2' '_f' num2str(SUBfac) '_Exp' expnum '.mat']);
        hh = waitbar(0,'Extract Epitoms...');
        for i = 1:size(mu,4)
            ordPatches = likelyPatches( mu(:,:,:,i), pcl(:,:,i), patchSize, Nep );
            histo_entropy = zeros(size(ordPatches,4),1);
            for p = 1:size(ordPatches,4)
                histo_entropy(p)  = entropy(ordPatches(:,:,1,p)) + entropy(ordPatches(:,:,2,p)) + entropy(ordPatches(:,:,3,p));
                % 		histo_entropy(p)  = entropy(rgb2gray(ordPatches(:,:,:,p)));
            end
            [dum ord] = sort(histo_entropy,'descend');
            ord = ord(1:Nth);
            epitexture(i).patch(:,:,:,1:Nth) = ordPatches(:,:,:,ord); % even
            for p = 1:Nth
                img_hsv     =   rgb2hsv(ordPatches(:,:,:,ord(p)));
                epitexture(i).hist(:,p)=   [hist(raster(img_hsv(:,:,1)),[0:1/NBINs(1):1]) ...
                    hist(raster(img_hsv(:,:,2)),[0:1/NBINs(2):1])...
                    hist(raster(img_hsv(:,:,3)),[0:1/NBINs(3):1])];
                % 				epitexture(i).lbph(:,p) = lbp(rgb2gray(ordPatches(:,:,:,ord(p))), 1,8,mapping, 'h');
            end
            waitbar(i/(size(mu,4)),hh);
        end
        close(hh);
    case 'ETHZ3'
        load(['MAT/Epitome_ETHZ3' '_f' num2str(SUBfac) '_Exp' expnum '.mat']);
        hh = waitbar(0,'Extract Epitoms...');
        for i = 1:size(mu,4)
            ordPatches = likelyPatches( mu(:,:,:,i), pcl(:,:,i), patchSize, Nep );
            histo_entropy = zeros(size(ordPatches,4),1);
            for p = 1:size(ordPatches,4)
                histo_entropy(p)  = entropy(ordPatches(:,:,1,p)) + entropy(ordPatches(:,:,2,p)) + entropy(ordPatches(:,:,3,p));
                % 		histo_entropy(p)  = entropy(rgb2gray(ordPatches(:,:,:,p)));
            end
            [dum ord] = sort(histo_entropy,'descend');
            ord = ord(1:Nth);
            epitexture(i).patch(:,:,:,1:Nth) = ordPatches(:,:,:,ord); % even
            for p = 1:Nth
                img_hsv     =   rgb2hsv(ordPatches(:,:,:,ord(p)));
                epitexture(i).hist(:,p)=   [hist(raster(img_hsv(:,:,1)),[0:1/NBINs(1):1]) ...
                    hist(raster(img_hsv(:,:,2)),[0:1/NBINs(2):1])...
                    hist(raster(img_hsv(:,:,3)),[0:1/NBINs(3):1])];
                % 				epitexture(i).lbph(:,p) = lbp(rgb2gray(ordPatches(:,:,:,ord(p))), 1,8,mapping, 'h');
            end
            waitbar(i/(size(mu,4)),hh);
        end
        close(hh);
% 	case 'ETHZ3'
% 		
% 		% load epitexture from each images
% 		load('MAT/epitomo_SEQ3_secondoTestA')
% 		muA = mu;
% 		pclA = pcl;
% 		
% 		load('MAT/epitomo_SEQ3_secondoTestB')
% 		muB = mu;
% 		pclB = pcl;
% 		
% 		clear mu mu2 mu3 pcl pcl2 pcl3 mask_seq3 x stat % cleaning
% 		
% 		% class A
% 		hh = waitbar(0,'Distances computation Epitexture...');
% 		for i = 1:size(muA,4)
% 			ordPatches = likelyPatches( muA(:,:,:,i), pclA(:,:,i), patchSize, Nep );
% 			histo_entropy = zeros(size(ordPatches,4),1);
% 			for p = 1:size(ordPatches,4)
% 				histo_entropy(p)  = entropy(ordPatches(:,:,1,p)) + entropy(ordPatches(:,:,2,p)) + entropy(ordPatches(:,:,3,p));
% 				% 		histo_entropy(p)  = entropy(rgb2gray(ordPatches(:,:,:,p)));
% 			end
% 			[dum ord] = sort(histo_entropy,'descend');
% 			ord = ord(1:Nth);
% 			epitexture(2*i-1).patch(:,:,:,1:Nth) = ordPatches(:,:,:,ord); % even
% 			for p = 1:Nth
% 				img_hsv     =   rgb2hsv(ordPatches(:,:,:,ord(p)));
% 				epitexture(2*i-1).hist(:,p)=   [hist(raster(img_hsv(:,:,1)),[0:1/NBINs(1):1]) ...
% 					hist(raster(img_hsv(:,:,2)),[0:1/NBINs(2):1])...
% 					hist(raster(img_hsv(:,:,3)),[0:1/NBINs(3):1])];
% % 				epitexture(2*i-1).lbph(:,p) = lbp(rgb2gray(ordPatches(:,:,:,ord(p))), 1,8,mapping, 'h');
% 			end
% 			waitbar(i/(size(muA,4)*2),hh);
% 		end
% 		
% 		% class B
% 		for i = 1:size(muB,4)
% 			ordPatches = likelyPatches( muB(:,:,:,i), pclA(:,:,i), patchSize, Nep );
% 			for p = 1:size(ordPatches,4)
% 				histo_entropy(p)  = entropy(ordPatches(:,:,1,p)) + entropy(ordPatches(:,:,2,p)) + entropy(ordPatches(:,:,3,p));
% 			end
% 			[dum ord] = sort(histo_entropy,'descend');
% 			ord = ord(1:Nth);
% 			epitexture(2*i).patch(:,:,:,1:Nth) = ordPatches(:,:,:,ord);   % odd
% 			for p = 1:Nth
% 				img_hsv     =   rgb2hsv(ordPatches(:,:,:,ord(p)));
% 				epitexture(2*i).hist(:,p)=   [hist(raster(img_hsv(:,:,1)),[0:1/NBINs(1):1]) ...
% 					hist(raster(img_hsv(:,:,2)),[0:1/NBINs(2):1])...
% 					hist(raster(img_hsv(:,:,3)),[0:1/NBINs(3):1])];
% % 				epitexture(2*i).lbph(:,p) = lbp(rgb2gray(ordPatches(:,:,:,ord(p))), 1,8,mapping, 'h');
% 			end
% 			waitbar(i/(size(muB,4)*2),hh);
% 		end
% 		close(hh);
% 		
		
	case 'VIPeR'
		
		% load epitexture from each images
		load(['MAT/Epitome_VIPeR' '_f' num2str(SUBfac) '_Exp' expnum '.mat']);
		
		clear mu2 mu3 pcl2 pcl3 mask_seq3 x stat % cleaning
		
		% class A
		hh = waitbar(0,'Distances computation Epitexture...');
		for i = 1:size(mu,4)
			ordPatches = likelyPatches( mu(:,:,:,i), pcl(:,:,i), patchSize, Nep );
			histo_entropy = zeros(size(ordPatches,4),1);
			for p = 1:size(ordPatches,4)
				histo_entropy(p)  = entropy(ordPatches(:,:,1,p)) + entropy(ordPatches(:,:,2,p)) + entropy(ordPatches(:,:,3,p));
				% 		histo_entropy(p)  = entropy(rgb2gray(ordPatches(:,:,:,p)));
			end
			[dum ord] = sort(histo_entropy,'descend');
			ord = ord(1:Nth);
			epitexture(i).patch(:,:,:,1:Nth) = ordPatches(:,:,:,ord); % even
			for p = 1:Nth
				img_hsv     =   rgb2hsv(ordPatches(:,:,:,ord(p)));
				epitexture(i).hist(:,p)=   [hist(raster(img_hsv(:,:,1)),[0:1/NBINs(1):1]) ...
					hist(raster(img_hsv(:,:,2)),[0:1/NBINs(2):1])...
					hist(raster(img_hsv(:,:,3)),[0:1/NBINs(3):1])];
% 				epitexture(i).lbph(:,p) = lbp(rgb2gray(ordPatches(:,:,:,ord(p))), 1,8,mapping, 'h');
			end
			waitbar(i/(size(mu,4)),hh);
		end
		close(hh);
        
    case 'iLIDS'
        
        % load epitexture from each images
        load(['MAT/Epitome_iLIDS' '_f' num2str(SUBfac) '_Exp' expnum '.mat']);
        
        clear mu2 mu3 pcl2 pcl3 mask_seq3 x stat % cleaning
        
        % class A
        hh = waitbar(0,'Load Epitexture...');
        for i = 1:size(mu,4)
            ordPatches = likelyPatches( mu(:,:,:,i), pcl(:,:,i), patchSize, Nep );
            histo_entropy = zeros(size(ordPatches,4),1);
            for p = 1:size(ordPatches,4)
                histo_entropy(p)  = entropy(ordPatches(:,:,1,p)) + entropy(ordPatches(:,:,2,p)) + entropy(ordPatches(:,:,3,p));
                % 		histo_entropy(p)  = entropy(rgb2gray(ordPatches(:,:,:,p)));
            end
            [dum ord] = sort(histo_entropy,'descend');
            ord = ord(1:Nth);
            epitexture(i).patch(:,:,:,1:Nth) = ordPatches(:,:,:,ord); % even
            for p = 1:Nth
                img_hsv     =   rgb2hsv(ordPatches(:,:,:,ord(p)));
                epitexture(i).hist(:,p)=   [hist(raster(img_hsv(:,:,1)),[0:1/NBINs(1):1]) ...
                    hist(raster(img_hsv(:,:,2)),[0:1/NBINs(2):1])...
                    hist(raster(img_hsv(:,:,3)),[0:1/NBINs(3):1])];
                % 				epitexture(i).lbph(:,p) = lbp(rgb2gray(ordPatches(:,:,:,ord(p))), 1,8,mapping, 'h');
            end
            waitbar(i/(size(mu,4)),hh);
        end
        close(hh);
        
end