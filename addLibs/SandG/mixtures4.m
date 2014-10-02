function [bestk,bestpp,bestmu,bestcov,dl,countf] = mixtures4(y,kmin,kmax,regularize,th,covoption,foptions_DATA)
%
% Syntax:
% [bestk,bestpp,bestmu,bestcov,dl,countf] = mixtures3(y,kmin,kmax,regularize,th,covoption)
%
% Inputs:
%
% "y" is the data; for n observations in d dimensions, y must have
% d lines and n columns. 
%
% "kmax" is the initial (maximum) number of mixture components
% "kmin" is the minimum number of components to be tested
%
% "covoption" controls the covarince matrix
% covoption = 0 means free covariances
% covoption = 1 means diagonal covariances
% covoption = 2 means a common covariance for all components
% covoption = 3 means a common diagonal covarince for all components
%
% "regularize" is a regularizing factor for covariance matrices; in very small samples, 
% it may be necessary to add some small quantity to the diagonal of the covariances
%
% "th" is a stopping threshold
% 
% Outputs:
%
% "bestk" is the selected number of components
% "bestpp" is the obtained vector of mixture probabilities
% "bestmu" contains the estimates of the means of the components
%          it has bestk columns by d lines
% "bestcov" contains the estimates of the covariances of the components
%           it is a three-dimensional (three indexes) array 
%           such that bestcov(:,:,1) is the d by d covariance of the first
%           component and bestcov(:,:,bestk) is the covariance of the last
%           component
%
% "dl" contains the successive values of the cost function
% "countf" is the total number of iterations performed
%
% Written by:
%   Mario Figueiredo
%   Instituto Superior Tecnico
%   1049-001 Lisboa
%   Portugal
%
%   Email: mtf@lx.it.pt
%
%   2000, 2001, 2002
% 
verb=0; % verbose mode; change to zero for silent mode
bins = 40; % number of bins for the univariate data histograms for visualization
dl = []; % vector to store the consecutive values of the cost function 
[dimens,npoints] = size(y);

switch covoption
   case 0
      npars = (dimens + dimens*(dimens+1)/2); 
      %this is for free covariance matrices
   case 1
      npars = 2*dimens;   
      %this is for diagonal covariance matrices
   case 2 
      npars = dimens;   
      %this is for a common covariance matrix 
      %independently of its structure)
   case 3 
     npars = dimens; 
   otherwise 
     %the default is to assume free covariances
     npars = (dimens + dimens*(dimens+1)/2);                   
end       
nparsover2 = npars / 2;

% we choose which axes to use in the plot, 
% in case of higher dimensional data (>2)
% Change this to have other axes being shown
axis1 = 1;
axis2 = 2;

% kmax is the initial number of mixture components
k = kmax;

% indic will contain the assignments of each data point to
% the mixture components, as result of the E-step
indic = zeros(k,npoints);


% ---------------
% Initialization
% ---------------
%% MARIO
% Initialization: we will initialize the means of the k components 
% with k randomly chosen data points. Randperm(n) is a MATLAB function
% that generates random permutations of the integers from 1 to n.
randindex = randperm(npoints);
randindex = randindex(1:k);
estmu = y(:,randindex);

% the initial estimates of the mixing probabilities are set to 1/k
estpp = (1/k)*ones(1,k);
%% FINE MARIO


%% MARCO 2/07/2005
% Initialization using k_means, with K choosen as K_max
% step = floor(256/kmax);
% incr = 0;
% % for i=1:kmax
% %     tmp           = incr + step;
% %     MY_centres(i) = tmp;
% %     incr          = tmp;
% % end
% randindex = randperm(npoints);
% randindex = randindex(1:k);
% estmu = y(:,randindex);
% [MY_centres_k,MY_OPTIONS, My_pi_tmp, MY_ERRLOG] = kmeans(estmu', y', foptions_DATA);
% estmu = MY_centres_k';
% estpp = sum(My_pi_tmp,1)./npoints;
%% FINE MARCO 2/07/2005


% here we compute the global covariance of the data
globcov = cov(y');

for i=1:k
   % the covariances are initialized to diagonal matrices proportional
   % to 1/10 of the mean variance along all the axes.
   % Of course, this can be changed
   estcov(:,:,i) = diag(ones(1,dimens)*max(diag(globcov/10)));
end


% % we plot the data and the initial estimates
% if (dimens == 1)
%    [hh,xx] = hist(y,bins);
% %    barplot = bar(xx,hh/npoints/(xx(2)-xx(1)));
% %    set(barplot,'EdgeColor',[1 1 1]);
% %    set(barplot,'FaceColor',[0.75 0.75 0.75]);
% %    set(gca,'FontName','Times'); set(gca,'FontSize',14);
% %    drawnow
% %    hold on
%    miny = min(y);
%    maxy = max(y);
% %    plotgrid = miny:(maxy-miny)/200:maxy;
% %    mix = zeros(size(plotgrid));
%    for comp=1:k
%        mix = mix + estpp(comp)*uninorm(plotgrid,estmu(comp),estcov(comp));
%        plot(plotgrid,estpp(comp)*uninorm(plotgrid,estmu(comp),estcov(comp)),'k')
%    end
%    plot(plotgrid,mix,'Color','red','LineWidth',3);
%    text(plotgrid(5),0.9*max(mix),sprintf('k=%d',k),...
%        'FontName','Times','FontSize',18,'FontWeight','bold');
%    drawnow
%    hold off
% else
%    plot(y(axis1,:),y(axis2,:),'.','Color',[0.5 0.5 0.5]);
%    hold on
%    axis equal
%    set(gca,'FontName','Times','FontSize',22);
%    placex = get(gca,'Xlim'); placey = get(gca,'Ylim');
%    text(placex(1)+0.2,placey(2)-0.2,sprintf('k=%d',k),...
%       'FontName','Times','FontSize',22);
%    for comp=1:k
%        elipsnorm(estmu([axis1,axis2],comp),...
%                  estcov([axis1,axis2],[axis1,axis2],comp),2)
%    end
%    drawnow
%    hold off
% end

% having the initial means, covariances, and probabilities, we can 
% initialize the indicator functions following the standard EM equation
% Notice that these are unnormalized values
for i=1:k
   semi_indic(i,:) = multinorm(y,estmu(:,i),estcov(:,:,i));
   indic(i,:) = semi_indic(i,:)*estpp(i);
end

% we can use the indic variables (unnormalized) to compute the
% loglikelihood and store it for later plotting its evolution
% we also compute and store the number of components
countf = 1;
loglike(countf) = sum(log(sum(realmin+indic)));
dlength = -loglike(countf) + (nparsover2*sum(log(estpp))) + ...
          (nparsover2 + 0.5)*k*log(npoints);   
dl(countf) = dlength;
kappas(countf) = k;

% the transitions vectors will store the iteration
% number at which components are killed.
% transitions1 stores the iterations at which components are 
% killed by the M-step, while transitions2 stores the iterations
% at which we force components to zero.
transitions1 = []; 
transitions2 = []; 

% minimum description length seen so far, and corresponding 
% parameter estimates
mindl = dl(countf);
bestpp = estpp;
bestmu = estmu;
bestcov = estcov;
bestk = k;


k_cont = 1;    % auxiliary variable for the outer loop
while(k_cont)  % the outer loop will take us down from kmax to kmin components
cont=1;        % auxiliary variable of the inner loop
while(cont)    % this inner loop is the component-wise EM algorithm with the
               % modified M-step that can kill components
   if verb==1
      % in verbose mode, we keep displaying the minimum of the 
      % mixing probability estimates to see how close we are 
      % to killing one component
      disp(sprintf('k = %2d,  minestpp = %0.5g', k, min(estpp)));
   end   
   
   % we begin at component 1
   comp = 1;
   % ...and can only go to the last component, k. 
   % Since k may change during the process, we can not use a for loop
   while comp <= k
            % we start with the M step
            % first, we compute a normalized indicator function
            clear indic
            for i=1:k
                indic(i,:) = semi_indic(i,:)*estpp(i);
            end
            normindic = indic./(realmin+kron(ones(k,1),sum(indic,1)));
            
            % now we perform the standard M-step for mean and covariance
            normalize = 1/sum(normindic(comp,:));
            aux = kron(normindic(comp,:),ones(dimens,1)).*y;
            estmu(:,comp)= normalize*sum(aux,2);
            if (covoption == 0)|(covoption == 2)
               estcov(:,:,comp) = normalize*(aux*y') - estmu(:,comp)*estmu(:,comp)' ...
                                  + regularize*eye(dimens);
            else
               estcov(:,:,comp) = normalize*diag(sum(aux.*y,2)) - diag(estmu(:,comp).^2) ;
            end
            if covoption == 2
               comcov = zeros(dimens,dimens);
               for comp2 = 1:k
                   comcov = comcov + estpp(comp2)*estcov(:,:,comp2);
               end
               for comp2 = 1:k
                   estcov(:,:,comp2) = comcov;
               end
            end
            if covoption == 3
               comcov = zeros(dimens,dimens);
               for comp2 = 1:k
                   comcov = comcov + estpp(comp2)*diag(diag(estcov(:,:,comp2)));
               end
               for comp2 = 1:k
                   estcov(:,:,comp2) = comcov;
               end
            end
                       
            
            % this is the special part of the M step that is able to
            % kill components
            estpp(comp) = max(sum(normindic(comp,:))-nparsover2,0)/npoints;
            estpp = estpp/sum(estpp);
            
            % this is an auxiliary variable that will be used the 
            % signal the killing of the current component being updated
            killed = 0;
            
            % we now have to do some book-keeping if the current component was killed
            % that is, we have to rearrange the vectors and matrices that store the
            % parameter estimates
            if estpp(comp)==0
                killed = 1;
                % we also register that at the current iteration a component was killed
                transitions1 = [transitions1 countf];
                if comp==1
                      estmu = estmu(:,2:k);
                      estcov = estcov(:,:,2:k);
                      estpp = estpp(2:k);
                      semi_indic = semi_indic(2:k,:);
                else 
                     if comp==k
                        estmu = estmu(:,1:k-1);
                        estcov = estcov(:,:,1:k-1);
                        estpp = estpp(1:k-1);
                        semi_indic = semi_indic(1:k-1,:);
                     else
                        estmu = [estmu(:,1:comp-1) estmu(:,comp+1:k)];
                        newcov = zeros(dimens,dimens,k-1);
                        for kk=1:comp-1
                            newcov(:,:,kk) = estcov(:,:,kk);
                        end
                        for kk=comp+1:k
                            newcov(:,:,kk-1) = estcov(:,:,kk);
                        end
                        estcov = newcov;
                        estpp = [estpp(1:comp-1) estpp(comp+1:k)];
                        semi_indic = semi_indic([1:comp-1,comp+1:k],:);
                     end
                end
            
                % since we've just killed a component, k must decrease
                k=k-1; 
             end    
             
             if (covoption == 2)|(covoption == 3)
                for kk = 1:k
                   semi_indic(kk,:) = multinorm(y,estmu(:,kk),estcov(:,:,kk)); 
                end
             end
             
             if killed==0
                % if the component was not killed, we update the corresponding
                % indicator variables...
                semi_indic(comp,:) = multinorm(y,estmu(:,comp),estcov(:,:,comp)); 
                % ...and go on to the next component
                comp = comp + 1;  
             end
             % if killed==1, it means the in the position "comp", we now
             % have a component that was not yet visited in this sweep, and
             % so all we have to do is go back to the M setp without 
             % increasing "comp"
             
    end % this is the end of the innermost "while comp <= k" loop 
        % which cycles through the components
        
    % increment the iterations counter            
    countf = countf + 1;
    if countf>100
%         fprintf('Bad data!\n');
        cont=0;
    end
    
    clear indic
    clear semi_indic
    for i=1:k
        semi_indic(i,:) = multinorm(y,estmu(:,i),estcov(:,:,i));
        indic(i,:) = semi_indic(i,:)*estpp(i);
    end

    if k~=1
       % if the number of surviving components is not just one, we compute
       % the loglikelihood from the unnormalized assignment variables
       loglike(countf) = sum(log(realmin+sum(indic)));
    else
       % if it is just one component, it is even simpler
       loglike(countf) = sum(log(realmin+indic));
    end
    
    % compute and store the description length and the current number of components
    dlength = -loglike(countf) + (nparsover2*sum(log(estpp))) + ...
              (nparsover2 + 0.5)*k*log(npoints); 
    dl(countf) = dlength;
    kappas(countf) = k;
    
    % compute the change in loglikelihood to check if we should stop
    deltlike = loglike(countf) - loglike(countf-1);
    if (verb~=0)
        disp(sprintf('deltaloglike/th = %0.7g', abs(deltlike/loglike(countf-1))/th));
    end      
    if (abs(deltlike/loglike(countf-1)) < th)
       % if the relative change in loglikelihood is below the threshold, we stop CEM2
       cont=0;
    end

  end % this end is of the inner loop: "while(cont)"
  
%   figure 
%   if dimens == 1
%    [hh,xx] = hist(y,bins);
%    barplot = bar(xx,hh/npoints/(xx(2)-xx(1)));
%    set(barplot,'EdgeColor',[1 1 1]);
%    set(barplot,'FaceColor',[0.75 0.75 0.75]);
%    set(gca,'FontName','Times'); set(gca,'FontSize',14);
%    drawnow
%    hold on
%    miny = min(y);
%    maxy = max(y);
%    plotgrid = miny:(maxy-miny)/200:maxy;
%    mix = zeros(size(plotgrid));
%    for comp=1:k
%        mix = mix + estpp(comp)*uninorm(plotgrid,estmu(comp),estcov(comp));
%        plot(plotgrid,estpp(comp)*uninorm(plotgrid,estmu(comp),estcov(comp)),'k')
%    end
%    plot(plotgrid,mix,'Color','red','LineWidth',3);
%    text(plotgrid(5),0.9*max(mix),sprintf('k=%d',k),...
%        'FontName','Times','FontSize',18,'FontWeight','bold');
%    drawnow
%    hold off
%  else
%   plot(y(axis1,:),y(axis2,:),'.','Color',[0.5 0.5 0.5]);
%   placex = get(gca,'Xlim'); placey = get(gca,'Ylim');
%   text(placex(1)+0.2,placey(2)-0.2,sprintf('k=%d',k),...
%        'FontName','Times','FontSize',22);
%   hold on
%   axis equal
%   set(gca,'FontName','Times','FontSize',22);
%   for comp=1:k
%       elipsnorm(estmu([axis1,axis2],comp),estcov([axis1,axis2],[axis1,axis2],comp),2)
%   end
%   drawnow
%   hold off
%  end
%   

    
  % now check if the latest description length is the best; 
  % if it is, we store its value and the corresponding estimates 
  if dl(countf) < mindl
       bestpp = estpp;
       bestmu = estmu;
       bestcov = estcov;
       bestk = k;
       mindl = dl(countf);
  end
  
  % at this point, we may try smaller mixtures by killing the 
  % component with the smallest mixing probability and then restarting CEM2,
  % as long as k is not yet at kmin
  if k>kmin
     [minp indminp] = min(estpp);      
     % what follows is the book-keeping associated with removing one component
     if indminp==1
        estmu = estmu(:,2:k);
        estcov = estcov(:,:,2:k);
        estpp = estpp(2:k);
     else 
        if indminp==k
           estmu = estmu(:,1:k-1);
           estcov = estcov(:,:,1:k-1);
           estpp = estpp(1:k-1);
        else
           estmu = [estmu(:,1:indminp-1) estmu(:,indminp+1:k)];
           newcov = zeros(dimens,dimens,k-1);
           for kk=1:indminp-1
               newcov(:,:,kk) = estcov(:,:,kk);
           end
           for kk=indminp+1:k
               newcov(:,:,kk-1) = estcov(:,:,kk);
           end
           estcov = newcov;
           estpp = [estpp(1:indminp-1) estpp(indminp+1:k)];
        end
     end
     k=k-1;
     
     % we renormalize the mixing probabilities after killing the component
     estpp = estpp/sum(estpp);
     
     % and register the fact that we have forced one component to zero
     transitions2 = [transitions2 countf];
          
     %increment the iterations counter
     countf = countf+1;
     
     % ...and compute the loglikelihhod function and the description length
     clear indic
     clear semi_indic
     for i=1:k
        semi_indic(i,:) = multinorm(y,estmu(:,i),estcov(:,:,i));
        indic(i,:) = semi_indic(i,:)*estpp(i);
     end
     if k~=1
        loglike(countf) = sum(log(realmin+sum(indic)));
     else
        loglike(countf) = sum(log(realmin+indic));
     end
     dl(countf) = -loglike(countf) + (nparsover2*sum(log(estpp))) + ...
              (nparsover2 + 0.5)*k*log(npoints); 

     
     kappas(countf) = k;
 
   else %this else corresponds to "if k > kmin"
        %of course, if k is not larger than kmin, we must stop      
        k_cont = 0;
   end
   
end % this is the end of the outer loop "while(k_cont)" 

lastpp = estpp;
lastmu = estmu;
lastcov = estcov;

% % finally, we plot the results
% figure 
%   if dimens == 1
%    [hh,xx] = hist(y,bins);
%    barplot = bar(xx,hh/npoints/(xx(2)-xx(1)));
%    set(barplot,'EdgeColor',[1 1 1]);
%    set(barplot,'FaceColor',[0.75 0.75 0.75]);
%    set(gca,'FontName','Times'); set(gca,'FontSize',14);
%    drawnow
%    hold on
%    miny = min(y);
%    maxy = max(y);
%    plotgrid = miny:(maxy-miny)/200:maxy;
%    mix = zeros(size(plotgrid));
%    for comp=1:bestk
%        mix = mix + bestpp(comp)*uninorm(plotgrid,bestmu(comp),bestcov(comp));
%        plot(plotgrid,bestpp(comp)*uninorm(plotgrid,bestmu(comp),bestcov(comp)),'k')
%    end
%    plot(plotgrid,mix,'Color','red','LineWidth',3);
%    text(plotgrid(5),0.9*max(mix),sprintf('k=%d',bestk),...
%        'FontName','Times','FontSize',18,'FontWeight','bold');
%    drawnow
%    hold off
% else
% plot(y(axis1,:),y(axis2,:),'.','Color',[0.5 0.5 0.5]);
% hold on
% axis equal
% set(gca,'FontName','Times','FontSize',22);
% placex = get(gca,'Xlim'); placey = get(gca,'Ylim');
% text(placex(1)+1,placey(2)-1,sprintf('k=%d',length(bestpp)),...
%      'FontName','Times','FontSize',22);
% for comp=1:length(bestpp)
%     elipsnorm(bestmu([axis1,axis2],comp),bestcov([axis1,axis2],[axis1,axis2],comp),2)
% end
% drawnow
% hold off
% end 
% thefig=gcf;
% 
% % now, we plot the evolution of the description length
% % and signal where components where killed
% figure
% plot(dl,'LineWidth',2)
% title('Description Length')
% ax = axis;
% for i=1:length(transitions1)
%     line([transitions1(i) transitions1(i)],[ax(3) ax(4)])
% end
% for i=1:length(transitions2)
%     line([transitions2(i) transitions2(i)],[ax(3) ax(4)],'LineStyle',':')
% end
% figure(thefig)

