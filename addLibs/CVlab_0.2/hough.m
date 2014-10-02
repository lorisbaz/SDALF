function res=hough(im,RHO_MAX,THETA_MAX)
%HOUGH Compute 2D Hough tansform
%
% USE:			res=hough(im,RHO_MAX,THETA_MAX)
%
% Name:			hough
%
% Version:		
%              v2.0
%
% Author: 		Dimitrios Ioannou
%			dimitris@cyra.com
%
%
% Date:			
%			v.1	08/23/95
%			v.1.1	03/13/96
%        v.2.0 04/29/99
%
% Arguments:
%			im: is the input,binary, image. If the
%			image is not binary pixels having
%			non-zero values are considered as feature points.
%			RHO_MAX: is an integer number specifying
%			the rho quantization.
%			THETA_MAX: is an integer number
%			specifying the theta quantization
%
% Purpose:
%			perform the Hough Transform of a square binary
%			image. Significantly faster than version 1.
%
% Dependencies:
%			None
%
% Example:		v=hough(im,256,256)
%			input is the image im, and the
%			quantization is d_rho=X/256 and d_theta=pi/256
%			if the size of the image is 256 by 256
%			d_rho=1.
%
%			v is the number of votes in the
%			parameter space. v(i,j) is the number 
%			of votes for the strip having distance from 
%			the center of the image equal to 
%			(i-RHO_MAX/2)*d_rho (d_rho=X/RHO_MAX, the 
%			image is X by X	pixels),and its normal has 
%			angle j*d_theta,(d_theta=pi/THETA_MAX)
%
%			for a 256 by 256 image, the center of the
%			image is the center of the pixel (128,128)
%			i=1 => rho=(i-1-128)*d_rho=-128*d_rho
%			i=256 => rho=(i-1-128)*d_rho=127*d_rho
%			this essentially means that:
%			'the image is not symmetric around its center'.
%
% BUGS FIXES:
%        does not crash when there is one/zero feature points
%

if (nargin~=3 | nargout~=1)
	fprintf(1,'Correct use: res=hough2(im,RHO_MAX,THETA_MAX).\n');
	error('Exiting...\n');
end

[X,Y]=size(im);
if X~=Y
	fprintf(1,'Input image is not square.\n');
	error('Exiting...\n');
elseif rem(X,2)==1
	fprintf(1,'Input image size has to be even in pixels.\n');
	error('Exiting...\n');
end

%hack to make it work when there is one feature point
if sum(im(:)) == 1
   [i,j]=find(im);
   
   is=i-1:i+1;
   js=j-1:j+1;
   
   is=is(find(is>0&is<size(im,1)));
   js=js(find(js>0&js<size(im,2)));
   
   im_new=im;
   im_new(is,js)=1;
   
   im2=im_new;
   im2(i,j)=0;
   
   h_new=hough3(im_new,RHO_MAX,THETA_MAX);
   h2=hough3(im2,RHO_MAX,THETA_MAX);
   
   res=h_new-h2;
   
   return;
elseif sum(im(:)) == 0
   
   res=zeros(RHO_MAX,THETA_MAX);
   
   return;
end

%%tic

d_rho=X/RHO_MAX;
d_theta=pi/THETA_MAX;

theta=0:d_theta:pi-d_theta;

smat=sin(theta);
cmat=cos(theta);

fprintf('Finding feature points.\n');
[x,y]=find(im);

fprintf('Translating so the origin is in the middle of the image.\n');
fprintf('Doing the Hough Transform.\n');
h1=((y-Y/2-1) * smat + (x-X/2-1) * cmat )/d_rho;

fprintf('Rounding.\n');
h3=round(h1+(RHO_MAX+1)/2);

clear h1 im x y


% HACK TO MAKE IT FASTER!
%
%new stuff from here

% efficient counting, instead of using a for loop as in hough1
%
%
% an improvement in terms of speed
%
% here are the steps:
%
% 1) each column of h3 array is sorted 
%
% h3 array contains the calculated rho(s), each column has constant
%      theta
%
% 2) expand the image temp and calculate if the differences are
% are different than 0
% 
% if they are it means a new rho, if not same rho
%
% 3) find nonzero values of difference
%
% from that we find the rhos (from K)
% and the thetas (from j) 
% and the votes (from i)
%
% 4) having all these, it's easy to built the Hough matrix
% 
%

% step 1
temp=flipud(sort(h3));

% step 2
difference=(diff([max(temp)+1;temp])~=0);

fprintf('calc sorting-differences time:\n');

% step 3
K=find(difference);

[i,j]=find(difference);
rho=temp(K);

votes=zeros(size(i));
for counter=1:length(j)-1
   
   if j(counter)==j(counter+1)
      votes(counter)=i(counter+1)-i(counter);
   else
      votes(counter)=size(h3,1)+1-i(counter);
   end
   
end
votes(length(j))=size(h3,1)+1-i(length(j));

% step 4
res=zeros(RHO_MAX,THETA_MAX);

for counter=1:length(j)
   
   if rho(counter)>-1 & rho(counter)<RHO_MAX
      res(rho(counter)+1,j(counter))=votes(counter);
   end
end
%%toc
