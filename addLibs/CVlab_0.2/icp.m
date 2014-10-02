function G = icp(model, data)
%ICP Iterative Closet Point algorithm.
%
%G = icp(model, data) dati 2 insiemi di dati model e data;
%entrambi corrispondenti ad un'unica forma ma espressi in 
%sistemi di riferimento diversi; viene calcolata la trasformazione
%rigida che minimizza la distanza tra i due insiemi.

% Author: A. Fusiello, 2003

eps = 0.00001; 
maxiter = 50;

G=eye(4);
res = inf;
prevres = 0;
i=0; % iterations counter
while ((abs(res-prevres)> eps) & i < maxiter)
    i = i+1;
    prevres = res;
    % apply current transformation bringing data onto model 
    dataREG = rigid(G(1:3,:),data);
    % compute closest points (in the model) and residuals
    [res,modelCP] = closestp(dataREG,model);
    % compute incremental tranformation
    GI = absolute(modelCP,dataREG);
    G = [GI; 0 0 0 1] * G;
end

G=G(1:3,:);


function  [res,modelcp] = closestp(data,model);
%CLOSESTP calcola i closest points.
%
%[res,modelcp] = closestp(data,model) calcola i closest point  dati due 
%insiemi di dati, e restituisce la distanza media tra i closest point.

data=data';
model=model';

modelcp = ones(size(data));
mindist = inf*ones(size(data,1),1);

for i = 1:size(data,1)
    for j = 1:size(model,1)    
        d = norm(model(j,:) - data(i,:));
        if d < mindist(i)
            mindist(i)=d;
            modelcp(i,:) =  model(j,:);
        end
    end
end

% apply X84
location = median(mindist);
scale = 5.2 * median(abs(mindist-location));
% set points to NaN in order to discard them
I = find(abs(mindist-location) > scale);
modelcp(I,:)=NaN;
% compute average distance for inliers
J = find(abs(mindist-location) <= scale);
res = mean(mindist(J));
modelcp=modelcp';
