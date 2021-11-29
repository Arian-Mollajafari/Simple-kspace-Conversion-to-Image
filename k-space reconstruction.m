clc;
clear;
clf;
load kspace;
whos;

[ydim, xdim] = size(kspace);
graymap = [0:1/255:1;0:1/255:1;0:1/255:1]';
subplot(1,2,1);
image(real(kspace),'CdataMapping','scaled');
axis square;
set(gcf,'Color','White');
colormap(graymap);
title('k-space real part');

subplot(1,2,2);
image(imag(kspace),'CdataMapping','scaled');
axis square;
title('k-space imaginary part');
%% Center of K-Space
Ymid = ydim/2;
Xmid = xdim/2;
kspace_center = kspace( Ymid-32:Ymid+31, Xmid-32:Xmid+31 );

subplot(1,2,1);
image(real(kspace_center),'CdataMapping','scaled');
axis square;
set(gcf,'Color','White');
colormap(graymap);
title('k-space real part');

subplot(1,2,2);
image(imag(kspace_center),'CdataMapping','scaled');
axis square;
title('k-space imaginary part');
%% The Image
clf;
ReconImage = ifftn(fftshift(kspace));

subplot(1,2,1);
image(abs(ReconImage), 'CdataMapping','direct');
axis image;
title('Reconstructed Image');
colormap(graymap);

% A lower resolution image
LR_image = ifftn(fftshift(kspace_center));
subplot(1,2,2);
image(abs(LR_image/64), 'CdataMapping','direct');
title('Low Resolution Image');
axis image;
%% Gaussian k-space 
% Create the gaussian
x=0:511;
SW = 32;
variance = (xdim/SW)^2;
g=repmat(exp(-(x-256).^2/variance),512,1);
Sen_map = g' .* g;

clf;
subplot(1,2,1);
image(abs(Sen_map),'CDataMapping','scaled');
title('Sensitivity Map');
axis image;
colormap(graymap);


Gaussian_kspace = Sen_map .* kspace;

Smoothed_image = ifftn(fftshift(Gaussian_kspace));

subplot(1,2,2);
image(abs(Smoothed_image), 'CdataMapping','scaled');
title('Smoothed Image');
axis image;

%% Higher variannce
SW = 2;
variance = (xdim/SW)^2;
g=repmat(exp(-(x-256).^2/variance),512,1);
Sen_map = g' .* g;

subplot(1,2,1);
image(abs(Sen_map),'CDataMapping','scaled');
title('Sensitivity Map');
axis image;
colormap(graymap);

Gaussian_kspace = Sen_map .* kspace;
Smoothed_image = ifftn(fftshift(Gaussian_kspace));

subplot(1,2,2);
image(abs(Smoothed_image), 'CdataMapping','scaled');
title('Less Smoothed Image');
axis image;

%% Construct high pass filter

SW = 2.5;
variance = (xdim/SW)^2;
g=repmat(exp(-(x-256).^2/variance),512,1);
Sen_map = g' .* g;

HP_map = fftshift(Sen_map);

subplot(1,2,1);
image(abs(HP_map),'CDataMapping','scaled');
title('Sensitivity Map');
axis image;
colormap(graymap);

HPkspace = HP_map .* kspace;
HighPassRecon = ifftn(fftshift(HPkspace));

subplot(1,2,2);
image(abs(HighPassRecon), 'CdataMapping','scaled');
title('High pass filtered image');
axis image;

%% Common artifacts

clear;
mygray = [0:1/255:1;0:1/255:1;0:1/255:1]';

load kspace;
kspaceC = kspace(129:384, 129:384);
clear kspace;
Spike = kspaceC;
Spike(120,120) = 1e7;
clf;

subplot(1,3,1);
imagesc(abs(ifftn(fftshift(kspaceC))));
colormap(mygray);
title('Original Image');
axis image;

subplot(1,3,2);
imagesc(abs(Spike));
title('kspace with Spike');
axis image;

subplot(1,3,3);
imagesc(abs(ifftn(fftshift(Spike))));
title('Image with spike in kspace');
axis image;

%% Two Spikes
subplot(1,3,1);
imagesc(abs(ifftn(fftshift(kspaceC))));
colormap(mygray);
title('Original Image');
axis image;

Spike(130,150) = 1e7;
subplot(1,3,2);
imagesc((abs(Spike)));
title('kspace with Spike');
axis image;

subplot(1,3,3);
imagesc(abs(ifftn(Spike)));
title('Image with spikes in kspace');
axis image;