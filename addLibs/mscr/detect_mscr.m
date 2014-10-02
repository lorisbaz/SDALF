%function [mvec,pvec[,arate,elist2]]=detect_mscr(image[,pars])
%
% Maximally Stable Colour Region detector
%
% IMAGE  Input image RGB or grey-scale double
% PARS   Optional parameter struct. Fields:
%          min_margin (default 0.0015)
%          timesteps  (default 200)
%          min_size   (default 60)
%          ainc (default 1.05)
%          filter_size, odd (default 3)
%          n8flag (default 0)
%          normfl (default 1)
%          blurfl (default 0)
%          verbosefl (default 1)
%
% MVEC   Moment vector list (6xN)
% PVEC   Property (colour) vector list (DxN)
% ARATE  Additional region properties (6xN)
%        Rows: area_rate,amin,d0,dn,t0,tn
% ELIST  Used edges
%        Rows: edge_value,x,y,dirn
%
% NOTE: DETECT_MSCR is a MEX function, and thus needs
%       to be compiled for your platform.
%
%Per-Erik Forssen, June 2007
function [mvec,pvec,arate,elist2]=detect_mscr(image,pars)

error(sprintf('Matlab cannot find a binary for this platform.\nHave you compiled detect_mscr yet?'));
