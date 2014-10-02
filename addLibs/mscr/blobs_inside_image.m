%function keep_list=blobs_inside_image(mvec,imsz)
%
% Remove blobs which are partially outside box imsz
%
% IMSZ  [rows cols] allowed region for blobs
%
% Returns an index list with indices of blobs inside image
%
%Per-Erik Forssen, Oct 2003

function [keep_list,BB]=blobs_inside_image(mvec,imsz)

[dummy,Nbl]=size(mvec);
keepfl=ones(1,Nbl);
COV=mvec([4 5 5 6],:);
M=mvec(2:3,:);

for k=1:Nbl,
  [E,D]=eig(sqrtm(reshape(COV(:,k),[2 2])));              
  D=2*D;                                % r=2*sqrt(lambda)
  % Find maximum extent in each direction
  xmax=sqrt((D(1,1)*E(1,1))^2+(D(2,2)*E(1,2))^2);
  ymax=sqrt((D(1,1)*E(2,1))^2+(D(2,2)*E(2,2))^2);
  % Blob inside image?
  keepfl(k)=((M(1,k)-xmax)>1)&((M(1,k)+xmax)<imsz(2))&...
  ((M(2,k)-ymax)>1)&((M(2,k)+ymax)<imsz(1));
  BB(:,:,k) = [M(:,k),[xmax;ymax]];
end

keep_list=find(keepfl);
