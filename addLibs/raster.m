function Is = raster(I)

switch length(size(I))
    case 2
        Is = I(:);
    case 3
        for i=1:3
            temp = I(:,:,i);
            Is(:,i) = temp(:);
        end
    otherwise
        fprintf('Unknown format\n');
        return
end
