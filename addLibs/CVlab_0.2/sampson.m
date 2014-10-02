function  d = sampson(F, m1, m2);
% SAMPSON compute Sampson error
%
% sampson(F, m1, m2)  valuate the first order approximation of the geometric error
% of the fit of F  with respect to a set of matched points m1, m2.
%
% Returns an approximation of the **squared** geometric distance from the
% 4d joint-space point [m1;m2] to the F manifold. 
%
% Author: Extracted and adapted from a piece of code by Peter Kovesi 

m1 = [m1;ones(1,size(m1,2))];
m2 = [m2;ones(1,size(m2,2))];

lng = size(m1,2);

m2tFm1 = zeros(1,lng);
for n = 1:lng
    m2tFm1(n) = m2(:,n)'*F*m1(:,n);
end

Fm1 = F*m1;
Ftm2 = F'*m2;

% Evaluate distances
d =  m2tFm1.^2 ./ (Fm1(1,:).^2 + Fm1(2,:).^2 + Ftm2(1,:).^2 + Ftm2(2,:).^2);


%------------------------------------------------------