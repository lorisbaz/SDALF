%% Re-identification 
% Loris Bazzani 
% Michela Farenzena

clear all;
close all;
clc; warning off;

% Additional libs loading
addpath(genpath('addLibs/'));


%% PARAMETRES
% Imgs param (for different tests)
SUBfac  = 1;    % subsampling factor <==== CHANGE TO 0.5 WHEN USING ETHZ DATASET!!!
H = 128*SUBfac; W = 64*SUBfac; % NORMALIZED dimensions
% symmetries var.
val    = 4;  
delta = [H/val W/val]; % border limit (in order to avoid the local minimums at the image border) 
varW    = W/5; % variance of gaussian kernel (torso-legs)
alpha = 0.5;
% HSV hist param 
NBINs   = [16,16,4]; % hsv quantization
% MSCR
parMSCR.min_margin	= 0.003; %  0.0015;  % Set margin parameter
parMSCR.ainc		= 1.05;
parMSCR.min_size	= 15;
parMSCR.filter_size	= 3;
parMSCR.verbosefl	= 0;  % Uncomment this line to suppress text output
% dynamic MSCR clustering parameters
kmin=1;
kmax=5;
regularize = 1;
th = 1e-2;
covoption=1;
% Matching param.
pyy = 0.4;
pcc = 0.6;
pee = 0.5;

% maximum number of images per signatures
MAXCLUSTER = 10; % <= THIS NUMBER OF IMAGES WILL BE SPLIT (N/2 for gallery, N/2 for probe)

% max texturized patch param.
N       = 30;	 % number of patches
fac     = [0,0;		% patches dim. (head is NOT USED!)
    12,10;			%(torso)
    10,8]*SUBfac;	%(legs)
var     = 3*SUBfac;     % dimension variance
NTrans      = 20;       % number of transformations (along y axis)
DIM_OP      = [35,35]*SUBfac;  % dim. local cross-correlation
Thresh_entr = 13;     % entropy thresholdings
Thresh_CC	= 0.4;		% crosscorrelation thresholding
mapping     = getmapping(8,'u2');			  % LBP param. (for CLUSTERING)
MSpar       = [10,25,45*ones(1,mapping.num)]; % MeanShift bandwidth

% head detection (ONLY IF USED)
% DIMW    = 24*SUBfac;  % variance of radial gaussian kernel (head)
% h_range = [0.0 0.1];  % skin det. param
% s_range = [0.3 0.6]; 


%% user choises (NO PARAMETERS!)
dataname    = 'iLIDS';	% dataset name 'iLIDS' 'ETHZ1' 'ETHZ2' 'ETHZ3'
expnum		= '012';	% exp number (for multiple experiments)
maskon		= 1;		% use or not the masks
dethead     = 0;		% enable/disable head detection
plotY       = 0;		% view partial results
maxplot		= 40;			% for visualization
loaddata    = 1;
basepath	= '/media/EA9C-FFBF/datasets/Tagging/'; % dataset path

switch dataname
	case 'iLIDS'
		avoid_inds = [];
%            11,14,33,42,60,67,68,91,92,96,127,128,147,148,149,165,195,210,...
% 			212,217,218,236,248,259,260,290,303,310,328,427,439,39,110,120,166,...
% 			167,178];  % not very good detected person  (I remove it!) % []; %
		dir_e   =   [basepath 'i-LIDS_Pedestrian'];
		dir_e_list = dir(strcat(dir_e,'/*.jpg'));
	case 'ETHZ1'
		avoid_inds = [];
		dir_e   =   [basepath 'dataset_ETHZ/seq1'];
		dir_e_list = dir(strcat(dir_e,'/*.png'));
		load('MAT/DynPedsModel_ETHZ1.mat') % ped(d).dynmodels contains the indexes of dyn model for the ped d
	case 'ETHZ2'
		avoid_inds = [];
		dir_e   =   [basepath 'dataset_ETHZ/seq2'];
		dir_e_list = dir(strcat(dir_e,'/*.png'));
		load('MAT/DynPedsModel_ETHZ2.mat') % ped(d).dynmodels contains the indexes of dyn model for the ped d
	case 'ETHZ3'
		avoid_inds = [];
		dir_e   =   [basepath 'dataset_ETHZ/seq3'];
		dir_e_list = dir(strcat(dir_e,'/*.png'));
		load('MAT/DynPedsModel_ETHZ3.mat') % ped(d).dynmodels contains the indexes of dyn model for the ped d

end
reg  	=   makecform('srgb2lab');
variance = var;

%% Init
n_img       =   size(dir_e_list,1);
permit_inds = setdiff(1:n_img,avoid_inds);
% permit_inds = 1:200; % ONLY for DEBUGGING!
n_img       =   length(permit_inds);
if plotY
    h1 = figure;
end
search_range_H  =   [delta(1),H-delta(1)];
search_range_W  =   [delta(2),W-delta(2)];

% Dataset inds extraction (selection of 2 views for each ped)
dset = SetDataset(permit_inds,dir_e_list);
if strcmp(dataname,'iLIDS')
	for i = 1:length(dset) % for each ped
		ped(i).dynmodels = dset(i).globalindex;
	end
end

if loaddata,
    
    %% Segmentatio FG/BG loading (+ morfological operations)
    if maskon
        load(['MAT/aleMasks_' dataname]) % load the ale masks
    end
    
    hwait = waitbar(0,'Mask processing...');
    for i=1:length(permit_inds)
        if maskon
            msk{permit_inds(i)} = imresize(msk{permit_inds(i)},[H,W]);
            mask_fin(:,:,i) = imfill(msk{permit_inds(i)},'holes');
        else
            mask_fin(:,:,i) = ones(H,W);
        end
        waitbar(i/length(permit_inds),hwait)
    end
    close(hwait)
    
    %% Dataset Loading
    fprintf([dataname ' Dataset LOADING...'])
    hwait = waitbar(0,'Dataset LOADING...');
    for i=1:length(permit_inds)
        img     =   imread(strcat(dir_e,'/',dir_e_list(permit_inds(i)).name));
        img     =   imresize(img,[H,W]);    % normalization 64x128
        dataset(:,:,:,i) = img;
        waitbar(i/length(permit_inds),hwait)
    end
    fprintf('OK \n');
    close(hwait)
end

%% Division in 3 part and kernel map computation
namefile = ['MAT/mapKern_div3_' dataname '_f' num2str(SUBfac) '_Exp' expnum '.mat'];
A = dir(namefile);
if ~isempty(A)
	load(namefile)
	fprintf('Division in 3 part and kernel map computation LOADING... OK\n')
else
	fprintf('Division in 3 part and kernel map computation COMPUTATION... ')
    tic
	mapkern_div3;
    tt(1) = toc;
	save(namefile,'MAP_KRNL','TLanti','BUsim','LEGsim','HDanti',...
		'head_det','head_det_flag');     %% saving
	fprintf('OK \n');
end


%% Features Extraction
    % 1) MSCR
	namefile = ['MAT/MSCR_' dataname '_f' num2str(SUBfac) '_Exp' expnum '.mat'];
	A = dir(namefile);
    if ~isempty(A)
        load(namefile)
        fprintf('MSCR LOADING... OK\n')
    else
        fprintf('MSCR COMPUTATION... ')
        tic
        ExtractMSCR; 
        tt(2) = toc;
        save(namefile,'Blobs');     %% saving
        fprintf('OK \n');
    end 
    
    % 2) Weighted HSV histogram
	namefile = ['MAT/wHSV_' dataname '_f' num2str(SUBfac) '_Exp' expnum '.mat'];
	A = dir(namefile);
    if ~isempty(A)
        load(namefile)
        fprintf('Weighted HSV hist LOADING... OK\n')
    else
        fprintf('Weighted HSV hist COMPUTATION... ')
        tic
        EstimatewHSV;
        tt(3) = toc;
        save(namefile,'whisto2');     %% saving
        fprintf('OK \n');
    end
    
%     % 3) Epitomo
% 	fprintf('Epitexture hist COMPUTATION... ')
% 	ExtractEpitexture;
% 	fprintf('OK \n');

%  3) Representive-textured patch extraction
	namefile = ['MAT/TxPatch_' dataname '_f' num2str(SUBfac) '_Exp' expnum '.mat'];
	A = dir(namefile);
    if ~isempty(A)
        load(namefile)
        fprintf('Epitextures extraction LOADING... OK\n')
    else
        fprintf('Epitextures extraction COMPUTATION \n')
        tic
        ExtractTxpatch;
        tt(4) = toc;
        save(namefile,'max_txpatch');     %% saving
        fprintf('OK \n');
    end

    
for t = 1:10
	% component of dynamic feature selection
	fprintf(['Test ' num2str(t) '...'])
    ped_back = ped;
	for i = 1:length(dset) % for each ped
        in = randperm(length(ped_back(i).dynmodels)); % random selection of the data 4 dynamic feature
        ped(i).dynmodels = ped_back(i).dynmodels(in(1:min(length(ped_back(i).dynmodels), MAXCLUSTER)));
        ped(i).rnddd = 1:length(ped(i).dynmodels);
	end
%% Matching
    % 1) Matching MSCR
    
	MSCRmatch_DynVSDyn_Cl;
    
	% 2) Matching wHSV
    wHSVmatch_DynVSDyn;
        
	% 3) Matching epitomo
    Epitextmatch_DynvsDyn;
    
%% TESTING: Evaluation of our Method
	
    crossvalidation_Dyn;
	fprintf('Ok \n')
    ped = ped_back;
end
fprintf('AUC: %d     normalized AUC: %f \n',mean(stats.AUC), mean(stats.nAUC))
figure, plot(mean(stats.CMC,1)/size(stats.CMC,2)*100, 'b.-', 'Linewidth',1.2); title('Cumulative Matching Characteristic (CMC)'), grid on, axis([1,25,0,100]);
figure,plot(M,mean(stats.SRR,1),'o-'),title('Synthetic Recognition Rate (SRR)'), grid on, axis([min(M),max(M),0,100])
save(['RES/stats_' dataname '_f' num2str(SUBfac) '_Exp' expnum '.mat'],'stats','ped')