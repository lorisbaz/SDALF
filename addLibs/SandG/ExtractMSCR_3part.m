%% part-based MSCR detection
hwait = waitbar(0,'Extract blobs from boby and legs');
for i=1:length(permit_inds)

    ph(i) = double(HDanti(i));
    pb(i) = double(TLanti(i));

    Ia = squeeze(dataset(:,:,:,i));

    % Illuminance normalization
    [Ia] = illuminant_normalization(Ia);
    
    % masking
    B = double(mask_fin(:,:,i));
    Ia = double(Ia).*cat(3,B,B,B); % mask application

    % Equalization
%     [H,S,V] = rgb2hsv(uint8(Ia));
%     Ve = histeq(V(B==1)); Veq = V; Veq(B == 1) = Ve;
%     Ia = cat(3, H,S,Veq);
%     Ia = hsv2rgb(Ia);

    % part-based MSCR computation + outliers elimination
    [mah, pah] = detection(Ia,B,[1:ph(i)],1,0);
    [mab, pab] = detection(Ia,B,[ph(i)+1:pb(i)],1,0);
    [mal, pal] = detection(Ia,B,[pb(i)+1:H],1,0);

    mab(3,:) = mab(3,:)+ph(i);
    mal(3,:) = mal(3,:)+pb(i);

    mvec = [mah mab mal];
    pvec = [pah pab pal];
	pvec = pvec/256;
    Blobs(i).mvec = mvec;
    Blobs(i).pvec = pvec;

	% part 1
	MSCR(i,1).mvec = mah;
	MSCR(i,1).pvec = pah;
	% part 2
	MSCR(i,2).mvec = mab;
	MSCR(i,2).pvec = pab;
	% part 3
	MSCR(i,3).mvec = mal;
	MSCR(i,3).pvec = pal;
	
    waitbar(i/length(permit_inds),hwait);

end
close(hwait);

