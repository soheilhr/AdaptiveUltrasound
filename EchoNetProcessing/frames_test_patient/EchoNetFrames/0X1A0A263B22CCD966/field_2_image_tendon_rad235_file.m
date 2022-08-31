clear all
% addpath("Field_II")

field_init(-1); %initialize Field II
c = 1540;  %(m/s) global speed of sound in medium
f0 = 11.579e6;  %(Hz) center freq. (...use to define xducer response)
% alpha = (0.5)*100/(1e6);  %(dB/m/Hz) attenuation of medium
fs = 100e6;  %(Hz) sampling freq. of system
set_field('c',c);  %(m/s)
% set_field('att',alpha);
set_field('fs',fs);  %(Hz)
% Set up transducer geometry
%Set up transducer geometry and display transducer geometry:
no_elements = 601;  %number of physical elements.
width = 46.5e-6;  %(m) width of elements in x-direction TODO try first with lambda/2, then with lambda
height = 6e-3;  %(m) width of elements in y-direction
kerf = 10e-6;  %(m) distance between elements (in x-direction)
no_sub_x = 1;  %number of sub-divisions of elements in x-direction
no_sub_y = 1;  %number of sub-divisions of elements in y-direction
%Rfocus = 3e-3;
focus = [0,0,3e-3];  %(m) Fixed focus for array (x,y,z). Vector with three elements. 
% set at 3 cm

Tx = xdc_linear_array(no_elements, width, height, kerf, no_sub_x, no_sub_y, focus);
Rx = xdc_linear_array(no_elements, width, height, kerf, no_sub_x, no_sub_y, focus);
% show_xdc(Tx)
% Define the transducer impulse response as a gaussian-modulated sinusoidal pulse with fractional -6dB bandwidth of 0.50.
fracBW = .70;  %fractional bandwidth of the xducer
% Fractional bandwidth = FWHM/fc -- center frequency
tc = gauspuls('cutoff', f0, fracBW, -6, -40);  %cutoff time at -40dB, fracBW @ -6dB
t = -tc:1/fs:tc;  %(s) time vector centered about t=0
impResp = gauspuls(t,f0,fracBW);
xdc_impulse(Tx,impResp);
figure; plot(t*1e6,impResp); grid on; xlabel('t (\musec)');
title('Transducer Impulse Response')
% Set the transducer excitation pulse to be a delta function (i.e., we want the impulse response).
excitation = 1;  %(impulse)
xdc_excitation(Tx,excitation);  % Th contains a "pointer" to the transducer object
% peak transmitted pressure fieldf from WkspEx3_beam
% Defining the transmit sequence
% In Field II, it is often easier to move the scatterers rather than moving the transducer. (The xdc_* functions always center the transducer at the origin.)
% In this example, we will define a linear scan format with a beamspacing of , i.e., with Nyquist lateral sampling. We will move the transducer from  to  and image several scatterers.
wl = c / f0;  % wavelength
D = no_elements * (width + kerf) - kerf;  % aperture width

Dy = height;
Dx = (width+kerf)*no_elements-kerf;
z = 3e-3; % focal depth

x_res = wl*z/Dx;
y_res = wl*z/Dy;
z_res = wl/2;
voxel_resolution = ['voxel resolution = ', num2str(x_res*1e3), 'mm x',num2str(y_res*1e3), 'mm x', num2str(z_res*1e3), 'mm'];
disp(voxel_resolution)


BeamSpacing = 1/2 * wl * focus(3) / D;  % beam spacing at half of the lateral resolution (f # = last term)
BeamOris = -10*wl:BeamSpacing:10*wl;  % define beam origins - at face of transducer
BeamOris = BeamOris - mean(BeamOris);  % make sure the origins are centered about x=0
nbeams = length(BeamOris);
% Define Scaterrer as a Tendon
% First, define a scatterer at 3cm depth (i.e., 1cm shallow to the focus).
% % % move the scatterer bc its hard to move the transducer in Field 2 --

% zposx = (-4e-3:1e-5:4e-3)'; 
% zposz = (2e-3:.725e-5:7.8e-3)';
% ScatLoc = [zposx, 0*zposx, zposz];
x_size = 3/1000; % width of phantom
y_size = 3/1000; % transverse width of phantom
z_size = 3/1000; % height of phantom
total_no_voxels = round((x_size/x_res)*(y_size/y_res)*(z_size/z_res));
for j = 1:15
    no_scatterers(j,1) = total_no_voxels*j;
end

Scatterers_phantom = round(no_scatterers);
% Define scatterers
index=12;
N=round(Scatterers_phantom(index));

% Scatterers positions
xp= (-x_size/2 + (x_size/2-(-x_size/2)).*rand(N,1));  %[m]
yp= (-y_size/2 + (y_size/2-(-y_size/2)).*rand(N,1));  %{m}
zp= (-z_size/2 + (z_size/2-(-z_size/2)).*rand(N,1)) + z;  %[m]

ScatLoc2 = [xp yp zp];
ScatAmp2 = randn(size(zp));

zposx = linspace(-3e-3,3e-3, 101)';
%zposz = linspace(2e-3, 4e-3, 101)';
zposz = 3e-3*ones(length(zposx),1);
%y = zposx*0;
%[xx, yy, zz] = ndgrid(zposx, y, zposz);
% xx, yy, and zz will all be 5 x 1 x arrays. 
ScatLoc = [zposx, zposx*0, zposz];
for i = 1:2
    ScatLoc = cat(1,ScatLoc, [zposx+1e-4*rand(1,1)-1e-4*rand(1,1), zposx*0, zposz+ 1e-4*rand(1,1) - 1e-4*rand(1,1)+i*3e-5]);
end

for i = 1:2
    ScatLoc = cat(1,ScatLoc, [zposx-1e-4*rand(1,1)+1e-4*rand(1,1), zposx*0, zposz+ 1e-4*rand(1,1) - 1e-4*rand(1,1)-i*3e-5]);
end


%ScatLoc = cat(1, ScatLoc, [zposx, zposx*0, zposz+1e-4]);
%ScatAmp = ones(size(zposx));
ScatAmp = ones(size(ScatLoc,1),1); % 1001 x 1 arrays
% Execute a pulse-echo scattering simulation using calc_scat:
% Pre-allocate rfdata and t0 arrays
rfdata_all = cell(nbeams,1);
t0_all = zeros(nbeams,1);
xdc_dynamic_focus(Tx, -1, 0, 0);
% Loop over beams
for b = 1:nbeams
    ScatLoc_b = ScatLoc - [BeamOris(b),0,0]; % get distance 
    [rfdata_all{b}, t0_all(b)] = calc_scat(Tx,Rx,ScatLoc_b,ScatAmp);
end


rfdata_all_2 = cell(nbeams,1);
t0_all_2 = zeros(nbeams,1);
% Loop over beams
BeamSpacing = 1/2 * wl * focus(3) / Dx;  % beam spacing at half of the lateral resolution
BeamOris2 = -x_size/2:BeamSpacing:x_size/2;  % define beam origins
BeamOris2 = BeamOris2 - mean(BeamOris2);  % make sure the origins are centered about x=0
nbeams = length(BeamOris2);

for b = 1:nbeams
    ScatLoc_b2 = ScatLoc2 - [BeamOris2(b),0,0]; % get distance 
    [rfdata_all_2{b}, t0_all_2(b)] = calc_scat(Tx,Rx,ScatLoc_b2,ScatAmp2);
end
% Now, the pulse-echo RF data from each simulation has a different start time (t0) and a different length. We need to align all of the traces. To do so, find the earliest t0, and also find the latest sample out of all the beams.
% Get the first and last sample index of each beam
startsamps = zeros(nbeams,1);
lastsamps = zeros(nbeams,1);
for b = 1:nbeams
    startsamps(b) = round(t0_all(b) * fs);
    lastsamps(b) = startsamps(b) + length(rfdata_all{b}) - 1;
end
% Now we know the total number of samples
firstsamp = min(startsamps);
nsamps = max(lastsamps) - firstsamp + 1;
rfdata = zeros(nsamps, nbeams);
% Populate rfdata with the samples in the right location
for b = 1:nbeams
    idx_start = startsamps(b) - firstsamp + 1;
    idx_end = lastsamps(b) - firstsamp + 1;
    rfdata(idx_start:idx_end,b) = rfdata_all{b};
end

% Get the first and last sample index of each beam
% For speckle
startsamps2 = zeros(nbeams,1);
lastsamps2 = zeros(nbeams,1);
for b = 1:nbeams
    startsamps2(b) = round(t0_all_2(b) * fs);
    lastsamps2(b) = startsamps2(b) + length(rfdata_all_2{b}) - 1;
end
% Now we know the total number of samples
firstsamp2 = min(startsamps2);
nsamps2 = max(lastsamps2) - firstsamp2 + 1;
rfdata2 = zeros(nsamps2, nbeams);
% Populate rfdata with the samples in the right location
for b = 1:nbeams
    idx_start2 = startsamps2(b) - firstsamp2 + 1;
    idx_end2 = lastsamps2(b) - firstsamp2 + 1;
    rfdata2(idx_start2:idx_end2,b) = rfdata_all_2{b};
end

rfdata3 = rfdata + rfdata2;
% We now have the RF data for all beams on the same time axis. Convert samples to time and get the time vector. Then, convert the time vector to depth.

% Image

toffset = -2*tc;  % Time offset due to sim start + time lag due to Tx and Rx impulse responses
time = (firstsamp:max(lastsamps)) / fs + toffset; %(s) time vector, accounting for offsets
depth = time * c/2+.1e-3;

% Finally, to make an image, we take several steps:
% Detect the envelope (magnitude of the Hilbert transform along the sample dimension).
% Normalize the image with respect to the maximum value.
% OPTIONAL: Resample the image to make the pixel size more isotropic.
% Log compress the image, i.e. 
% Display with a desired dynamic range, say, 50 dB.

bimg = abs(hilbert(rfdata3));  % Detect the envelope
bimg = bimg / max(bimg(:));  % Normalize the image
upsampfactor = 8;
bimg = interp2(bimg, 1:1/upsampfactor:nbeams, (1:nsamps)');  % Optional: interpolate bimg laterally
bimg = db(bimg);  % Log-compress the image

imgx = BeamOris*1e3;  % x-positions
imgx = interp1(imgx, 1:1/upsampfactor:nbeams);  % Also have to upsample the x-positions to match the image
imgz = depth*1e3;
imagesc(imgx, imgz, bimg, [-50 0]); grid on; xlabel('Lateral (mm)'); ylabel('Depth (mm)'); axis image; colormap gray
ylim([2.5 3.5]);
