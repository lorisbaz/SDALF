function [S1,S2] = ExtractDatasetOrd(dset)

% extract train/probe data

S1 = [];
S2 = [];
for i=1:length(dset)
    if length(dset(i).index) > 1
        %in = randperm(length(dset(i).index));
        S1 = [S1; dset(i).index(1)];
        S2 = [S2; dset(i).index(2)];
    end
end

