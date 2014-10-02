function lines = flin1(I,NLINES,meth)
%FLIN1 Find lines in an image with Hough transform.
%   If NLINES>0 the best NLINES lines are returned; the local maxima
%   suppression method is precise but slow. 
%   If NLINES=0 (default) all the lines with more than MAXFRAC relative
%   consensus are returned (MAXFRAC = 1 returns all the lines); the local
%   maxima suppression method is fast, based on morphologhical operations,
%   but it's not very precise.
%   meth is passed to the EDGE function (see) and specifies the method to be used
%   for edge extraction. Default: meth = 'sobel'.
%   Output format is compatible with DRAW_LINES (see)
%
%   See also: HOUGH, EDGE, DRAW_LINES

% Author: Andrea Fusiello

%% Trova linee con Hough di Dimitrios Ioannou (dal sito Mathworks).
%% Bugs: Funziona solo con immagini quadrate, e la funzione chiamata hough3 non esiste.

if nargin == 1
    NLINES = 0;
    meth='sobel';
elseif nargin ==2
    meth='sobel';
end

tic

DILATEFRAC=.01; % amount of dilation
MAXFRAC=0.6; % soglia per i picchi della HT, come frazione del massimo
THETA_MAX = 360; % quantizzazione dell'angolo
RHO_MAX = 2*ceil(norm(size(I)-floor((size(I)-1)/2)-1))+3; % formula per RHO_MAX come Radon
WS = 5; %finestra per la soppressione dei massini locali


E = edge(I,meth);
figure, imshow(E)%, title ('Edge Image')  

%Compute the Hough transform of the edge image.
R = hough(E,RHO_MAX,THETA_MAX);
figure, imagesc(R); colormap(hot);
%title('Hough transform');


if NLINES > 0
    
    %%Find peaks in the linespace: method 1  
    % -----------------------------------------
    
    % sopprime massimi locali
    
    % add a very small random noise to break ties (TAPULLO)
    R = R + max(R(:))*1e-10*rand(size(R)) ;
    % pad R
    R = [R(:,end-WS+1:end), R,  R(:,1:WS)];
    A=colfilt(R,[WS WS],'sliding',@max);
    R = R(:,WS+1:end-WS);
    A = A(:,WS+1:end-WS);
    local_maxima = find(R<A);
    R(local_maxima) = min(R(:));
    
    
    % trova i picchi sopra soglia
    m = max(R(:));
    i = find(R>MAXFRAC*m);  
    %%[y,x] = find(R>MAXFRAC*m)  
    
    %Sort the output and pick the top lines
    [foo,ind] = sort(-R(i));
    k = i(ind(1:NLINES));  
    
    %Convert linear index into coordinates of the peaks.
    [y,x] = ind2sub(size(R),k);  
    
    % -----------------------------------------
elseif NLINES==0 
    
    %%Find peaks in the linespace: method 2
    % -----------------------------------------
    
    m = max(max(R));
    TH=im2bw(R/m,MAXFRAC);
    
    %Morphological dilation
    H1=dilate(TH,ones(round(DILATEFRAC*size(R))),1);
    
    %Labeling
    L=bwlabel(H1,8);
    n=max(L(:));
    
    %for the n lines found above we collect the indices into the HT matrix
    x=[]; y=[];
    for k=1:n,
        [r,c]=find(L==k);
        y=[y; mean(r)];
        x=[x; mean(c)];
    end
    
    % Nota: questa tecnica non funziona bene quando vi sono rette perfettamente
    % orizzontali; in quel caso ci sono picchi alle estremita' destra e sinistra della
    % trasformata e bisognerebbe che le dilatazione e le componenti connesse fossero 
    % eseguite con il bordo destro e sinistro coincidenti.
    
    % -----------------------------------------
end

figure, imagesc(R), colormap(gray);
%%numplot([x y]), 
%%title('Location of peaks');  

%Find the theta and rho values for the peak coordinates.  
t = (x-1)*pi/THETA_MAX;
r = (y-1-RHO_MAX/2)*size(I,1)/RHO_MAX;

%The line parameters have the coefficents of the equation Ax + By + C = 0
lines = [sin(t) cos(t) -r];  

%Transform the line from the center of the image to the upper left.  
cx = floor((size(I,2)+1)/2);
cy = floor((size(I,1)+1)/2);
lines(:,3) = lines(:,3) - lines(:,1)*cx - lines(:,2)*cy;  

toc
