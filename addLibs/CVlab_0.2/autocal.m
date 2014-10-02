function K = autocal(imagePoints)
%AUTOCAL Computes intrinsic parameters with autocalibration
%
% Algorithm: Medonca & Cipolla 1999
%
% Author: A. Fusiello, 2005

siz = size(imagePoints);

numberOfPoints = siz(2);
numberOfViews= siz(3);



% calcolo le n*(n-1)/2  matrici fondamentali di tutte le coppie di viste

for i=1:numberOfViews
  for j=i+1:numberOfViews
    
    % calcola la matrice fondamentale ij 
    Fs(:,:,i,j) = ...
	fm(imagePoints(:,:,i),imagePoints(:,:,j)); 
    
    % check fundamental matrix by computing residuals of m_j * F_ij * m_i = 0  
     err = trace(abs([imagePoints(:,:,j); ones(1,numberOfPoints)]' *  Fs(:,:,i,j) * ...
	 [imagePoints(:,:,i); ones(1,numberOfPoints)]))/numberOfPoints;
    %disp(sprintf('Fundamental matrix %d %d residual: %0.5g',i,j,err));
    Ws(i,j) = 1/err;
  end
  
end

% normalizza i pesi
Ws = Ws./sum(sum(Ws));


% parameters initial guess:  dim = [width,height]
dim = max(max(imagePoints,[],3),[],2);

a = [dim(1)+dim(2), dim(1)+dim(2), dim(1)/2, dim(2)/2];

% find intrinsic parameters by minimizing residuals with Nelder-Meads

a = fminsearch(@(x) hf_cost(x,Fs,Ws),a);

% a = fminunc(@(x) hf_cost(x,Fs,Ws),a);

% this is the resulting Intrinsic parameters mateix 

K = [
  a(1)	0     a(3)
  0	a(2)  a(4)
  0	0      1];



return



function cost = hf_cost(params,Fs,Ws)
% Huang anf Faugeras objective function for autocalibration

% F(:,:,i,j) is the 3x3 fundamental matrix between views i and j
% params is the a 5 vector containing the intrinsic parameters
% of all the cameras

% Author: Andrea Fusiello


siz = size(Fs);
n = siz(3);


cost = 0.0;

Aj = [ params(1)  0  params(3)
      0 params(2) params(4)
      0 0 1];
    
Ai = [ params(1)  0 params(3)
      0 params(2) params(4)
      0 0 1];
    
for i=1:n
  for j=i+1:n 
    E = Aj'* Fs(:,:,i,j) * Ai;
    D = E*E';
    cost =  cost +  Ws(i,j) * (2.0 * trace(D*D') - trace(D)^2 )/ trace(D)^2;
  end
end

return




function cost= mc_cost(params,Fs,Ws)
% Medonca & Cipolla cost function for autocalibration

% F(:,:,i,j) is the 3x3 fundamental matrix between views i and j
% params is the a 5 vector containing the intrinsic parameters
% of all the cameras

% Author: Andrea Fusiello


siz = size(Fs);


n = siz(3);

Aj = [ params(1)  0  params(3)
         0 params(2) params(4)
         0 0 1];

Ai = [ params(1)  0 params(3)
        0 params(2) params(4)
        0 0 1];

cost = 0.0;

for i=1:n
  for j=i+1:n
    
    
    E = Aj'* Fs(:,:,i,j) * Ai;
    
    [U,S,V] = svd(Aj'* Fs(:,:,i,j) * Ai);
    
    cost = cost + Ws(i,j) * (S(1,1)-S(2,2))/(S(1,1)+S(2,2));
    
  end
end

return


