%function [Bimg] = draw_blobs(mvec,pvec,rows,cols[,bkgr])
%
% Visualise detected blob features
%
% MVEC   Moment vector list (6xN)
% PVEC   Property (colour) vector list (DxN)
% ROWS   Desired number of rows in output bitmap
% COLS   Desired number of columns in output bitmap
% BKGR   Background colour vector (default green, i.e. [0 1 0])
%
% BIMG   Blob image (approximating ellipses for regions)
%
% NOTE: DRAW_BLOBS is a MEX function, and thus needs
%       to be compiled for your platform.
%
%Per-Erik Forssen, June 2007
function Bimg = draw_blobs(mvec,pvec,rows,cols,bkgr)

error(sprintf('Matlab cannot find a binary for this platform.\nHave you compiled draw_blobs yet?'));
