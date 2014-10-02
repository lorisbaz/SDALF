function [In] = illuminant_normalization(IMG)

gray = [98.0 98.0 98.0 ];

I = double(IMG);
R = I(:,:,1);
G = I(:,:,2);
B = I(:,:,3);

Rn = R./mean(R(:))*gray(1);
Gn = G./mean(G(:))*gray(2);
Bn = B./mean(B(:))*gray(3);

In = Rn;
In(:,:,2) = Gn;
In(:,:,3) = Bn;
