%% ******** sequent 21 compound transmitting, 21 receiving ********

clear;

%% load data
load('0.mat')
data = pt;

%% ******** data process ********
%% parameter settings
N_elements = 64;
pitch = 0.4/1000;
Ntim    = 2560;
cs      = 1540;
f0      = 3e6;
fs      = 62.5e6/2;
xele    = (0:N_elements-1) * pitch;
XLIM    = xele(end);
YLIM    = 0.045;
densx   = 4;
dx      = pitch/densx;
dy      = cs/(20*f0);
x       = 0 : dx : XLIM;
y       = 0.005 : dy : YLIM;
Theta = (-10:10)*1/180*pi;

%% synthetic aperture beamforming, using Delay and Sum
RF = int32(zeros(length(y),length(x)));
temp = RF;
for i = 1:64
    temp = DAScompound64(data(1:Ntim*(length(Theta)),2*i-1),single(xele),single(x),single(y),int8(i-1));
    RF = RF + temp;
end

%% linear interpolation
dy2   = cs/(40*f0);
y2    = 0.005 : dy2 : YLIM;
midrf = RF;
clear RF
RF2 = int32(zeros(length(y2),size(midrf,2)));

for m = 1:size(RF2,2)
    RF2(:,m) = interp1(y,double(midrf(:,m)),y2,'linear');
end

%% ******** plot ********
Bmode = 20*log10(abs( hilbert(RF2(:,:) )) +1 );
Bmode = int8(Bmode);
figure; imagesc(x,y2,Bmode);
axis equal
axis([0 0.0252 0.005 0.045]) ;
colormap gray
caxis([70, 110]);