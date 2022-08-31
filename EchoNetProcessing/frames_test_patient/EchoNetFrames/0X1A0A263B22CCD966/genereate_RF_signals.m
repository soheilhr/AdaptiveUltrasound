% Beamform signals from a phased array
param = getparam('P4-2v');

width = pi/3; % width angle in rad
txdel = txdelay(param,0,width); % in s

load('test_phantom_patient1_bd10.mat')

param.fs = 4*param.fc; % sampling frequency

RF = simus(sca_z, sca_z, ampl, txdel, param);

% Demodulate the RF signals
IQ = rf2iq(RF,param);
% create a 256 x 256 80 degree wide polar grid with IMPOLGRID
[x,z] = impolgrid([256 256],10e-2,pi/3,param);

% Beamform the I/Q signals: 
IQb = das(IQ,x,z,txdel,param);

% create the b mode image
B = bmode(IQb,30); % log-compressed B-mode image with a -30 dB range


% display the b mode image
pcolor(x*100,z*100,B)
c = colorbar;
c.YTick = [0 255];
c.YTickLabel = {'-30 dB','0 dB'};
colormap gray
title('A simulated ultrasound image')
ylabel('[cm]')
shading interp
axis equal ij tight
set(gca,'XColor','none','box','off')
