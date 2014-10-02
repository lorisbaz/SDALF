% VIPER

cam_a = '../VIPeR/cam_a/';
cam_b = '../VIPeR/cam_b/';
cam_a_list = dir(strcat(cam_a,'/*.bmp'));
cam_b_list = dir(strcat(cam_b,'/*.bmp'));

load('RES/stats_VIPeR_f1_Exp003.mat');

fid = fopen('RES/VIPeR_partitions.txt', 'w');
for i=1:10
    
    fprintf(fid, '%% Experiment Number %d \n\n', i);
    
    S = stats.index(i,:);
    S = sort(S, 'ascend');
    
    for j=1:length(S)
    
        fprintf(fid,'%d - %s  %s \n', j, ['cam_a/' cam_a_list(S(j)).name], ['cam_b/' cam_b_list(S(j)).name]);
        
    end
        
    VIPeRpartitions(i).cam_a = {cam_a_list(S).name};
    VIPeRpartitions(i).cam_b = {cam_b_list(S).name};
    
    fprintf(fid,'\n \n');
    
end

fclose(fid);





%% iLIDS

cam_e = '../datasets/i-LIDS_Pedestrian/';
cam_e_list = dir(strcat(cam_e,'/*.jpg'));


load('RES/stats_iLIDS_f1_Exp003.mat');
fid = fopen('RES/iLIDS_gallerysets.txt', 'w');
for i=1:10
     fprintf(fid, '%% Experiment Number %d \n\n', i);
    S = stats.S1(i,:);
    S = sort(S, 'ascend');
    
     for j=1:length(S)
    
        fprintf(fid,'%d - %s \n', j, [cam_e_list(S(j)).name]);
        
    end
        fprintf(fid,'\n \n');
   
    iLIDSparitions(i).models = {cam_e_list(S).name};
  
end
fclose(fid);




%% ETHZ1

cam_e = '../datasets/dataset_ETHZ/seq1/';
cam_e_list = dir(strcat(cam_e,'/*.png'));


load('RES/stats_ETHZ1_f0.5_Exp003.mat');
fid = fopen('RES/ETHZ1_gallerysets.txt', 'w');

for i=1:10
     fprintf(fid, '%% Experiment Number %d \n\n', i);
    S = stats.A1(i,:);
    S = sort(S, 'ascend');
     for j=1:length(S)
    
        fprintf(fid,'%d - %s \n', j, [cam_e_list(S(j)).name]);
        
    end
        fprintf(fid,'\n \n');
   
    ETHZ1paritions(i).models = {cam_e_list(S).name};
  
end
fclose(fid);

%% ETHZ2

cam_e = '../datasets/dataset_ETHZ/seq2/';
cam_e_list = dir(strcat(cam_e,'/*.png'));


load('RES/stats_ETHZ2_f0.5_Exp003.mat');
fid = fopen('RES/ETHZ2_gallerysets.txt', 'w');

for i=1:10
    fprintf(fid, '%% Experiment Number %d \n\n', i);
    S = stats.A1(i,:);
    S = sort(S, 'ascend');
    ETHZ2paritions(i).models = {cam_e_list(S).name};
    for j=1:length(S)
    
        fprintf(fid,'%d - %s \n', j, [cam_e_list(S(j)).name]);
        
    end
        fprintf(fid,'\n \n');
end

fclose(fid);

%% ETHZ3

cam_e = '../datasets/dataset_ETHZ/seq3/';
cam_e_list = dir(strcat(cam_e,'/*.png'));


load('RES/stats_ETHZ3_f0.5_Exp003.mat');
fid = fopen('RES/ETHZ3_gallerysets.txt', 'w');
for i=1:10
    fprintf(fid, '%% Experiment Number %d \n\n', i);
    S = stats.A1(i,:);
    S = sort(S, 'ascend');
    ETHZ3paritions(i).models = {cam_e_list(S).name};
   for j=1:length(S)
    
        fprintf(fid,'%d - %s \n', j, [cam_e_list(S(j)).name]);
        
    end
        fprintf(fid,'\n \n');
end
fclose(fid);
