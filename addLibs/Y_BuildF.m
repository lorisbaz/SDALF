function [F] = Y_BuildF(data,coords,d)
% Y_BuildF Feature extractor.
%   F = Y_BuildF(data) returns the F values extracted on the pedestrian data
%   as described in [1]

%   References:
%      [1]O. Tuzel, F. Porikli, P. Meer: Pedestrian Detection via
%      Classification on Riemannian Manifolds, Pattern Analysis and Machine
%      Intelligence, IEEE Transactions on : Accepted for future publication
%          

%   Copyright 2008 Marco Cristani 
%   $Revision: 1.0 $  $Date: 2008/07/14$
% PS: must be converted in a function before to use into Z_PWatch 

%% MAIN
%---------------------- 
[size_arg] = size(data);
if length(size_arg)<2
    error('Data stored has to be W*H*num_images\n');
    return
end
W           =   size_arg(1);
H           =   size_arg(2);

[x,y]     =   meshgrid((coords(1):coords(3))',(coords(2):coords(4))'); 

img               =   double(data(:,:)) ;
F(:,:,1)          =   double(x)';
F(:,:,2)          =   double(y)';
[gh,gv]           =   gradient(img);
d_img_h           =   abs(gh);
d_img_v           =   abs(gv);
F(:,:,3)          =   d_img_h; %horizontal
F(:,:,4)          =   d_img_v; %vertical
F(:,:,5)          =   sqrt(gh.^2+gv.^2);
[gh]              =   gradient(gh);
F(:,:,6)          =   abs(gh);
[gh,gv]           =   gradient(gv);
F(:,:,7)          =   abs(gv);
F(:,:,8)          =   atan2(abs(d_img_h),(abs(d_img_v)));

%--- test uguaglianza metodi di calcolo
% F1          =   atan2(abs(d_img_h),(abs(d_img_v)));
% F2          =   atan(abs(d_img_h)./(abs(d_img_v)));
% A = F2;
% A = A(:);
% A(d_img_h == 0) = 0;
% A(d_img_v == 0 & d_img_h == 0) = 0;%A(d_img_v == 0 & d_img_h == 0) = pi/2; E` LA GESTIONE PASSATA
% A = reshape(A,W,H);
% figure;imagesc(A == F1)    
    