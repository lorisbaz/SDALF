function dset = SetDataset(permit_inds,dir_e_list)
% dset = SetDataset(permit_inds,dir_e_list)
%
% extract the index of each view of a single pedestrian
% dset is indexed by the pedestrian ID and contains the relative and global 
% indexes of the different views (imgs into the dataset)

j = 0; idpersonold = 0;

for i=1:length(permit_inds)
    
    n = dir_e_list(permit_inds(i)).name;
    idperson = str2num(n(1:4));
    nimage = str2num(n(5:7));
    
    if (idperson ~= idpersonold)
        j = j+1;
        dset(j).globalindex = permit_inds(i);
        dset(j).index = i;
        idpersonold = idperson;
    else
        dset(j).globalindex = [dset(j).globalindex permit_inds(i)];
        dset(j).index = [dset(j).index i];
    end
    
end