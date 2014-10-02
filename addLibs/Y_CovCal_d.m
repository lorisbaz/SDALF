function C = Y_CovCal_d(pxy,Qxy,R,d,R_main)
% C = Y_CovCal(pxy,Qxy,R) returns covariance descriptors C
% Y_CovarianceCalculation Reads tag of INRIA dataset.
% C = Y_CovCal_d(pxy,Qxy,R,d)
%
% In:
% pxy := first order tensor of integral images WxHxd
% Qxy := second order tensor of integral images WxHxd
% R := 4-tuple of the region (in absolute or relative coordinate)
% Out:
% C := covariance descriptors
%
%   References:
%      [1]O. Tuzel, F. Porikli, P. Meer: Pedestrian Detection via
%      Classification on Riemannian Manifolds, Pattern Analysis and Machine
%      Intelligence, IEEE Transactions on : Accepted for future publication
%          

%   Copyright 2009 (Marco Cristani), (Michela Farenzena), Diego Tosato
%   $Revision: 1.0 $  $Date: 2009/04/14$

% --- test
%S    = R_main(3)*R_main(4); % normalization term
%% Covariance matrix  computation for R 
  if(R(1) == 1 && R(2) == 1) % caso particolare, sono in alto a sx; 
    S    = R(3)*R(4); % normalization term
    to_transp = reshape(pxy(R(3),R(4),:),d,1);
    C = 1/(S-1).*(reshape(Qxy(R(3),R(4),:,:),d,d) -...
        1/S.*(to_transp*to_transp'));
  elseif R(1) == 1
      S   = (R(3) - R(1) + 1)*(R(4) - R(2) + 1);% normalization term  
      first_fact = reshape(Qxy(R(3),R(4),:,:) - Qxy(R(3),R(2)-1,:,:),d,d);
      to_transp = reshape((pxy(R(3),R(4),:) -  pxy(R(3),R(2)-1,:)),d,1);
      C = (1/(S-1)).*(first_fact-(1/S).*(to_transp*to_transp'));
  elseif R(2) == 1
      S    = (R(3) - R(1) + 1)*(R(4) - R(2) + 1);% normalization term
      first_fact = reshape(Qxy(R(3),R(4),:,:) -  Qxy(R(1)-1,R(4),:,:),d,d);
      to_transp = reshape((pxy(R(3),R(4),:) - pxy(R(1)-1,R(4),:)),d,1);
      C = (1/(S-1)).*(first_fact-(1/S).*(to_transp*to_transp'));
  else
      S    = (R(3) - R(1) + 1)*(R(4) - R(2) + 1);% normalization term
      first_fact = reshape(Qxy(R(3),R(4),:,:) + Qxy(R(1)-1,R(2)-1,:,:) - Qxy(R(3),R(2)-1,:,:) -  Qxy(R(1)-1,R(4),:,:),d,d);
      to_transp = reshape((pxy(R(3),R(4),:) + pxy(R(1)-1,R(2)-1,:) - pxy(R(1)-1,R(4),:) - pxy(R(3),R(2)-1,:)),d,1);
      C = (1/(S-1)).*(first_fact-(1/S).*(to_transp*to_transp'));
  end

  