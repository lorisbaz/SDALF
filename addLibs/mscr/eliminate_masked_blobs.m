%function keep_list=blobs_inside_image(mvec,imsz, Mask)
%
% Remove blobs composed by more than 70% of background pixels
%
% IMSZ  [rows cols] allowed region for blobs
%
% Returns an index list with indices of blobs inside image
%
%Per-Erik Forssen, Oct 2003

function [mvec, pvec] =eliminate_masked_blobs(mvec,pvec, imsz, Mask)

[dummy,Nbl]=size(mvec);
keepfl=ones(1,Nbl);
COV=mvec([4 5 5 6],:);
M=mvec(2:3,:);

keep_list = [];

for k=1:Nbl,
  %fprintf('%d \n', k);  
  [E,D]=eig(sqrtm(reshape(COV(:,k),[2 2])));              
  %D=2*D;                                % r=2*sqrt(lambda)
  % Find maximum extent in each direction
  xmax=2*sqrt((D(1,1)*E(1,1))^2+(D(2,2)*E(1,2))^2);
  ymax=2*sqrt((D(1,1)*E(2,1))^2+(D(2,2)*E(2,2))^2);
  
  
  xl = ceil((M(1,k)-xmax)-1); if xl<=0, xl= 1; end
  xh = floor((M(1,k)+xmax)-1); if xh>imsz(2), xh=imsz(2); end
  
  yl = ceil ((M(2,k)-ymax)-1); if yl<=0, yl= 1; end
  yh = floor((M(2,k)+ymax)-1); if yh>imsz(1), yh=imsz(1); end
  
  A = 0.25*E * inv(D./2) *E';
  
%   EI = zeros(imsz);
%   for x=xl:xh
%       for y=yl:yh
%           EI(y,x) = ([x-M(1,k) y-M(2,k)] * A * [x-M(1,k) y-M(2,k)]' < 1); 
%       end
%   end
%   
%   figure(1), subplot(1,2,1), imagesc(EI);
%   Bimg = draw_blobs(mvec(:,k),pvec(:,k),imsz(1),imsz(2),[0 255 0]');
%   subplot(1,2,2), image(uint8(Bimg));
%   pause
  
  XY = [repmat(xl:xh, 1, (yh-yl)+1); vec(repmat(yl:yh, (xh-xl)+1, 1))'];
  Mk = repmat(M(:,k), 1, size(XY,2));
  
  InsideEllipse = diag(single((XY-Mk)' * A * (XY-Mk) < 4));
  InsideEllipse = reshape(InsideEllipse,  (xh-xl)+1, (yh-yl)+1)';
  
  Areak = length(find(InsideEllipse));
  %AreaMasked = length(intersect(find(InsideEllipse), find(Mask(yl:yh,xl:xh)<200)));
  
  %if (AreaMasked/Areak) < 0.5 && (mvec(1,k) < 2*imsz(1)*imsz(2)/3)
  %  keep_list = [keep_list k];
  %end
  
end

% if isempty(keep_list)
%     m = find(mvec(1,:) == min(mvec(1,:)));
%     keep_list = m;
% end
% 
% mvec = mvec(:,keep_list);
% pvec = pvec(:,keep_list);
% 
%  Bimg = draw_blobs(mvec,pvec,imsz(1),imsz(2),[0 1 0]');
%  figure, image(Bimg), axis image
% % 
