%% dataset renamer
ETHZpath = '/Users/bazzani/Pictures/dataset_ETHZ/';

dirfolds = dir([ETHZpath '/seq*']);
for i = 1:length(dirfolds) % for each seq.
	dirpeds = dir([ETHZpath dirfolds(i).name '/p*']);
	
	for j = 1:length(dirpeds) % for each ped
		dirview = dir([ETHZpath dirfolds(i).name '/' dirpeds(j).name '/frame*']);
		if j<10
			suffixp = '000';
		elseif j<100
			suffixp = '00';
		elseif j<1000
			suffixp = '0';
		else
			suffixp = [];
		end
			
		for v = 1:length(dirview) % for each view
			if v<10
				suffixv = '00';
			elseif v<100
				suffixv = '0';
			else
				suffixv = [];
			end
		
			I = imread([ETHZpath dirfolds(i).name '/' dirpeds(j).name '/' dirview(v).name]);
% 			imagesc(I),pause
			
			imwrite(I,[ETHZpath 'new' num2str(i) '/' suffixp num2str(j) suffixv num2str(v) '.png']);
		end
	end
end