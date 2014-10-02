%% Re-identification
% Loris Bazzani 
% Michela Farenzena

clear all;
close all;
clc; warning off;

% Additional libs loading
addpath(genpath('addLibs/'));


%% PARAMETRES
% Imgs parameters (for different tests)
SUBfac  = 1;    % subsampling factor <==== CHANGE TO 0.5 WHEN USING ETHZ DATASET!!!
H = ceil(128*SUBfac); W = ceil(64*SUBfac); % NORMALIZED dimensions

% symmetries parameters
val    = 4;  
delta = [H/val W/val]; % border limit (in order to avoid the local minimums at the image border) 
varW    = W/5; % variance of gaussian kernel (torso-legs)
alpha = 0.5;

% HSV hist parameters 
NBINs   = [16,16,4]; % hsv quantization

% MSCR parameters
parMSCR.min_margin	= 0.003; %  0.0015;  % Set margin parameter
parMSCR.ainc		= 1.05;
parMSCR.min_size	= 15;
parMSCR.filter_size	= 3;
parMSCR.verbosefl	= 0;  % Uncomment this line to suppress text output

% Matching parameters
pyy = 0.4;
pcc = 0.6;
pee = 0.5;

% RHSP parameters
N       = 30;	 % number of patches
fac     = [0,0;		% patches dim. (head is NOT USED!)
    12,10;			%(torso)
    10,8]*SUBfac;	%(legs)
variance     = 3*SUBfac;     % dimension variance
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
dataname    = 'iLIDS';	% dataset name 'iLIDS' 'VIPeR' 'ETHZ1' 'ETHZ2' 'ETHZ3'
expnum		= '007';	% exp number (for multiple experiments)
maskon		= 1;		% use or not the masks
dethead     = 0;		% enable/disable head detection
plotY       = 0;		% view partial results
maxplot		= 40;			% for visualization
basepath	= '/media/EA9C-FFBF/datasets/Tagging/'; % dataset path
loaddata    = 1;        % if you need to load the dataset

% Extract image names
switch dataname
	case 'iLIDS'
		avoid_inds = [];
		dir_e   =   [basepath 'i-LIDS_Pedestrian'];
		dir_e_list = dir(strcat(dir_e,'/*.jpg'));
	case 'VIPeR'
        W = 48*SUBfac; 
        varW    = W/5;
		avoid_inds = [];
		dir_e   =   [basepath 'VIPeRa'];
		dir_e_list = dir(strcat(dir_e,'/*.bmp'));
        load('cvpridx.mat');
	case 'ETHZ1'
		avoid_inds = [];
		dir_e   =   [basepath 'dataset_ETHZ/seq1'];
		dir_e_list = dir(strcat(dir_e,'/*.png'));
		selected_models = [];
	case 'ETHZ2'
		avoid_inds = [];
		dir_e   =   [basepath 'dataset_ETHZ/seq2'];
		dir_e_list = dir(strcat(dir_e,'/*.png'));
		selected_models = [];
	case 'ETHZ3'
		avoid_inds = [];
		dir_e   =   [basepath 'dataset_ETHZ/seq3'];
		dir_e_list = dir(strcat(dir_e,'/*.png'));
		selected_models = [];
end
reg  	=   makecform('srgb2lab');


%% Init
n_img       =   size(dir_e_list,1);
permit_inds = setdiff(1:n_img,avoid_inds);
n_img       =   length(permit_inds);
if plotY
    h1 = figure;
end
search_range_H  =   [delta(1),H-delta(1)];
search_range_W  =   [delta(2),W-delta(2)];

% Dataset inds extraction (selection of 2 views for each ped)
dset = SetDataset(permit_inds,dir_e_list);

if loaddata
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

%% Division in 3 parts and kernel map computation
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
    
     %   3) RHSP
	namefile = ['MAT/txpatch_' dataname '_f' num2str(SUBfac) '_Exp' expnum '.mat'];
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
   

   
   
%% Matching
    % 1) Matching MSCR
	namefile = ['MAT/MSCRmatch_' dataname '_f' num2str(SUBfac) '_Exp' expnum '.mat'];
	A = dir(namefile);
	if ~isempty(A)
		load(namefile)
		fprintf('MSCR matching LOADING... OK\n')
	else
		fprintf('MSCR matching COMPUTATION... ')
        tic
		MSCRmatch;
        tt(5) = toc;
		save(namefile,'final_dist_y','final_dist_color','final_mscr_dist')
		fprintf('OK \n');
	end
	
	% 2) Matching wHSV
	namefile = ['MAT/wHSVmatch_' dataname '_f' num2str(SUBfac) '_Exp' expnum '.mat'];
	A = dir(namefile);
	if ~isempty(A)
		load(namefile)
		fprintf('wHSV matching LOADING... OK\n')
	else
		fprintf('wHSV matching COMPUTATION... ')
        tic
		wHSVmatch;
        tt(6) = toc;
		save(namefile,'final_dist_hist')
		fprintf('OK \n');
	end
	
	% 3) Matching RHSP
	namefile = ['MAT/txpatchmatch_' dataname '_f' num2str(SUBfac) '_Exp' expnum '.mat'];
	A = dir(namefile);
    if ~isempty(A)
        load(namefile)
        fprintf('Epitextures matching LOADING... OK\n')
    else
        fprintf('Epitextures matching COMPUTATION... ')
        tic
        CompareEpitext;
        tt(7) = toc;
        save(namefile,'dist_epitext')
        fprintf('OK \n');
    end
        
    
%% TESTING: Evaluation of our Method
	crossvalidation;
	save(['RES/stats_' dataname '_f' num2str(SUBfac) '_Exp' expnum '.mat'],'stats')	