%% ******** sequent 5 transmitting,20 receiving***************

clear

%% load data
load('0.mat');
data = pt;

%% data preprocess

%% parameter settings
N_elements = 64;
c = 1540; 
pitch = 0.4/1000;
Ntim    = 2560;
cs      = c;
f0      = 3e6;
fs      = 62.5e6/2;
xele    = (0:N_elements-1) * pitch;
XLIM    = xele(end);
YLIM    = 0.1;
densx = 4;
dx      = pitch/densx;
dy      = cs/(20*f0);
x       = 0 : dx : XLIM;
y       = (0 : dy : YLIM );
Theta = (-10:10)*1/180*pi;

%% synthetic aperture beamforming, using Delay and Sum

RF = zeros(length(y),length(x));
RF(:,:) = DAScompound64(double(data(1:Ntim*(length(Theta)),1:2:128)),fs,cs,xele,x,y,Theta,N_elements);

%% cubic spline interpolation
dy2   = cs/(40*f0);
y2    = 0 : dy2 : YLIM;
midrf = RF;
clear RF
RF2 = zeros(length(y2),size(midrf,2));

for m = 1:size(RF2,2)
    RF2(:,m) = interp1(y,double(squeeze(midrf(:,m))),y2,'Spline');
end

Bmode = 20*log10(abs( hilbert(squeeze(RF2(:,:)) )) +1 );
figure; mesh(x,y2,Bmode)
axis equal
colormap gray
view([0 0 1])
% caxis([ prctile(Bmode, 5,'all') , prctile(Bmode, 99,'all') ])