function patches = ptcs_yTrans(patch,NTrans,wplot)
% Generate NTrans patches (affine) rotated wrt the y axis.

if nargin<3
    wplot = false;
end
if wplot
    figure, subplot(2,NTrans/2+1,1), imshow(patch)
end

% declarations
patches = cell(NTrans,1);
% theta = (rand(NTrans,1)*2-1)*60;
theta = linspace(-pi,pi,NTrans);

% centering the points to the simmetry axis (y)
x = -(size(patch,2)-1)/2;  y = -(size(patch,1)-1)/2;
M = [x,y,0;         % UP-LEFT point
    x,y+size(patch,1),0; % DOWN-LEFT point
    x+size(patch,2),y,0; % UP-RIGHT point
    x+size(patch,2),y+size(patch,1),0]'; % DOWN-RIGHT point

for t = 1:NTrans
    % rotation of the (corner) points
    Ry = [cos(theta(t)),0,sin(theta(t)); 0,1,0; -sin(theta(t)),0,cos(theta(t))];
    for i = 1:size(M,2)
        m(:,i) = Ry*M(:,i);
    end
    % cancel z components
    Mi = M(1:2,:);
    mi = m(1:2,:);
    % traslate the point into image coordinate
    Mi(1,:) = Mi(1,:) + size(patch,2)/2;  
    Mi(2,:) = Mi(2,:) + size(patch,1)/2;
    mi(1,:) = mi(1,:) + (size(patch,2)-1)/2;
    mi(2,:) = mi(2,:) + (size(patch,1)-1)/2;

    % computing the homography
    H = hm(Mi,mi);   
    
    % warping of the image by the homography
    [patches{t},bb,alpha] = imwarp(patch,H,'linear','valid');

    % cutting of the black border
    bbox = regionprops(bwlabel(alpha),'BoundingBox');
    if ~isempty(bbox)
        bbox = ceil(bbox.BoundingBox);
        patches{t} = patches{t}(bbox(2):bbox(2)+bbox(4)-1,bbox(1):bbox(1)+bbox(3)-1);

        if wplot
            subplot(2,NTrans/2+1,t+1),imshow(patches{t}), axis image, axis on, colormap gray
            title(num2str(theta(t)))
        end
    else
        t = t-1;
    end
end
