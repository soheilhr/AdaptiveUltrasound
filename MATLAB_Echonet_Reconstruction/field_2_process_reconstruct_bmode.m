clear all;
close all;
% addpath("Field_II") 
tic

field_init(-1); %initialize Field II
c = 1540;  %(m/s) global speed of sound in medium
f0 =5e6;  %(Hz) center freq. (...use to define xducer response)
alpha = (0.5)*100/(1e6);  %(dB/m/Hz) attenuation of medium
fs = 100e6;
%set_sampling(fs);  %(Hz) sampling freq. of system CHECK WITH DONGWOON
set_field('c',c);  %(m/s)
set_field('att',alpha);
set_field('fs',fs);  %(Hz)

% Set up transducer geometry
%Set up transducer geometry and display transducer geometry:
no_elements = 64;  %number of physical elements.
width = 320e-6;  %(m) width of elements in x-direction TODO try first with lambda/2, then with lambda
height = 6e-3;  %(m) width of elements in y-direction
kerf = 20e-6;  %(m) distance between elements (in x-direction)
no_sub_x = 1;  %number of sub-divisions of elements in x-direction INCREASE
no_sub_y = 1;  %number of sub-divisions of elements in y-direction
%Rfocus = 3e-3;
focus = [0,0,10e-2];  %(m) Fixed focus for array (x,y,z). Vector with three elements. 

Tx = xdc_linear_array(no_elements, width, height, kerf, no_sub_x, no_sub_y, focus);
Rx = xdc_linear_array(no_elements, width, height, kerf, no_sub_x, no_sub_y, focus); % Define Tx and Rx separately
%show_xdc(Tx)

% Define the transducer impulse response as a gaussian-modulated sinusoidal pulse with fractional -6dB bandwidth of 0.50.
fracBW = .70;  %fractional bandwidth of the xducer ; Fractional bandwidth = FWHM/fc -- center frequency
tc = gauspuls('cutoff', f0, fracBW, -6, -40);  %cutoff time at -40dB, fracBW @ -6dB
t = -tc:1/fs:tc;  %(s) time vector centered about t=0
impResp = gauspuls(t,f0,fracBW);
xdc_impulse(Tx,impResp);
xdc_impulse(Rx, impResp); % Setting impulse response for both Tx and Rx
% figure; plot(t*1e6,impResp); grid on; xlabel('t (\musec)');
% title('Transducer Impulse Response')


% Set the transducer excitation pulse to be a delta function (i.e., we want the impulse response).
excitation = 1;  %(impulse) NUNMBER OF CYCLES
xdc_excitation(Tx,excitation);  % Th contains a "pointer" to the transducer object. Set excitation for Tx only. 
% peak transmitted pressure fieldf from WkspEx3_beam


% Defining the transmit sequence
% In Field II, it is often easier to move the scatterers rather than moving the transducer. (The xdc_* functions always center the transducer at the origin.)
% In this example, we will define a linear scan format with a beamspacing of , i.e., with Nyquist lateral sampling. We will move the transducer 
% from x to y and image several scatterers.
wl = c / f0;  % wavelength
D = no_elements * (width + kerf) - kerf;  % aperture width

Dy = height;
Dx = (width+kerf)*no_elements-kerf;
z = 10e-2; % focal depth

x_res = wl*z/Dx;
y_res = wl*z/Dy;
z_res = wl/2;
voxel_resolution = ['voxel resolution = ', num2str(x_res*1e3), 'mm x',num2str(y_res*1e3), 'mm x', num2str(z_res*1e3), 'mm'];
disp(voxel_resolution)


% BeamSpacing = 1/2 * wl * focus(3) / D;  % beam spacing at half of the lateral resolution (f # = last term)
%BeamOris = -10*wl:BeamSpacing:10*wl;  % define beam origins - at face of transducer
%BeamOris = BeamOris - mean(BeamOris);  % make sure the origins are centered about x=0
% Use the entire transducer as the beam origin

% for nbeams, do one beam per degree. Calculated as 88-89 degrees, 
nbeams = 89; % NEED TO MAKE THIS NOT HARD CODED; GET SLOPES OF EACH AND CALCULATE TAN(THETA) = (M1 - M2)/(1Mm1*M2)

% Beam Spacing is based on theta -- the angle of the image
sector = 89 * pi/180;
d_theta = sector/nbeams;

theta = -sector/2;

% Scatterer is mat file from MUST code
load('test_phantom_patient1_bd5_new_axis.mat')
%[xx, yy, zz] = ndgrid(zposx, y, zposz);
%sca_x = linspace(-60e-3,60e-3, 20)';
%sca_z = linspace(20e-3, 120e-3, 20)';
y = 0*sca_x;
ScatLoc = [sca_x, y, sca_z];
%ScatLoc = [0,0,4e-2];
ScatAmp = ampl;
%ScatAmp = ones(size(ScatLoc,1),1); % 1001 x 1 arrays
%ScatAmp = 1;

% Execute a pulse-echo scattering simulation using calc_scat:

% Pre-allocate rfdata and t0 arrays
rfdata_all = cell(nbeams,1);
t0_all = zeros(nbeams,1);
image_data = zeros(800, nbeams);
%xdc_focus(Tx, -1, ScatLoc); % test ; later, loop through for each transmit focus 
%xdc_dynamic_focus(Rx, -1, 0, 0)

% xdc_focus_times(Tx,-1*ones(17535,1),zeros(17535,64));
% Loop over beams
for b = 1:nbeams
    %ScatLoc_b = ScatLoc - [0,0,0]; % get distance ; beam origin for all points is (0,0,0)
    focal_points = [ScatLoc(:,1)*sin(theta) ScatLoc(:,2) ScatLoc(:,3)]; 
    time = -1*ones(length(ScatLoc), 1);
    xdc_focus(Tx, -time, focal_points);
    xdc_dynamic_focus(Rx, -1, 0, 0);
    % Calculate the received response
    [rfdata_all{b}, t0_all(b)] = calc_scat(Tx,Rx, ScatLoc,ScatAmp);
    image_data(1:max(size(rfdata_all{b})), b) = rfdata_all{b}';
    %times(b) = t0_all{b};
    theta = theta + d_theta
end
figure
plot(image_data)

% Now, the pulse-echo RF data from each simulation has a different start time (t0) and a different length. We need to align all of the traces. 
% To do so, find the earliest t0, and also find the latest sample out of all the beams.
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

bimg = abs(hilbert(rfdata));  % Detect the envelope
bimg = bimg / max(bimg(:));  % Normalize the image
upsampfactor = 8;
bimg = interp2(bimg, 1:1/upsampfactor:nbeams, (1:nsamps)');  % Optional: interpolate bimg laterally
bimg = db(bimg);  % Log-compress the image
x_axis = linspace(-.075, .075, 112);
imgx = x_axis*1e3;  % x-positions
imgx = interp1(imgx, 1:1/upsampfactor:nbeams);  % Also have to upsample the x-positions to match the image
imgz = depth*1e3;
imagesc(imgx, imgz, bimg, [-50 0]); grid on; xlabel('Lateral (mm)'); ylabel('Depth (mm)'); axis image; colormap gray; 
toc
t = toc;
