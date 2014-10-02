function dset = removePeds(dset,Nv)
% Remove the peds with less than Nv views in order to sudivide the dataset
% in two class of dynamic features

j = 1;
ok = false;
for i = 1:length(dset)
	if length(dset(i).index) >= Nv
		dset_red(j) = dset(i);
		j = j + 1;
		ok = true;
	end	
end
if ok
	dset = dset_red;
end