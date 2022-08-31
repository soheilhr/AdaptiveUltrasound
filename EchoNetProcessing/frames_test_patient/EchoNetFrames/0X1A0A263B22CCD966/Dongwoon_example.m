% Execute a Field II simulation
% Field II (https://field-ii.dk/) is a program for simulating ultrasound imaging, written by Prof. Jørgen Arendt Jensen of the Technical University of Denmark. This MATLAB Live script gives a quick tutorial for beamforming a simple Field II point target simulation.
% Simulation parameters
% Define the simulation parameters.
clear all
addpath ../matBF/  % Add path to rtbf_mex
fs = 80e6;  % Sampling frequency [Hz]
c0 = 1540;  % Sound speed [m/s]
f0 = 8e6;   % Modulation frequency [Hz]
wl = c0 / f0;  % Wavelength [m]
txfn = 1;  % Transmit f-number (z/D)
rxfn = 1;  % Receive f-number (z/D)
% Set the field of view (FOV) extents
limx = [-5e-3 5e-3];
limy = [0 0];
limz = [12e-3 28e-3];
% Transducer definition
% We will simulate an L12-3v transducer with 128 elements, 8 MHz center frequency, and 80% fractional bandwidth.
% Define element geometry
xdc.nx = 128;                               % Number of elements in x
xdc.ny = 1;                                 % Number of elements in y
xdc.n = xdc.nx * xdc.ny;                    % Total number of elements
xdc.pitch_x = 200e-6;                       % Element pitch in x
xdc.pitch_y = 7.5e-3;                       % Element pitch in y
xdc.kerf_x = 25e-6;                         % Kerf in azim
xdc.kerf_y = 0;                             % Kerf in y
xdc.width_x = xdc.pitch_x - xdc.kerf_x;     % Element width in x
xdc.width_y = xdc.pitch_y - xdc.kerf_y;     % Element width in y
xdc.Dx = xdc.nx * xdc.pitch_x - xdc.kerf_x; % Aperture width in x
xdc.Dy = xdc.ny * xdc.pitch_y - xdc.kerf_y; % Aperture width in y
xdc.radius = inf;                           % No curvature
xdc.elevfoc = 20e-3;                        % Mechanical lens focus
% Compute element positions
xdc.xc = (1:xdc.nx)' * xdc.pitch_x;
xdc.yc = (1:xdc.ny)  * xdc.pitch_y;
xdc.xc = xdc.xc - mean(xdc.xc);  % Center at origin
xdc.yc = xdc.yc - mean(xdc.yc);  % Center at origin
xdc.xc = repmat(xdc.xc, 1, xdc.ny);
xdc.yc = repmat(xdc.yc, xdc.nx, 1);
xdc.zc = zeros(xdc.nx, xdc.ny);
xdc.pos = [xdc.xc(:)'; xdc.yc(:)'; xdc.zc(:)'];  % [3, nelems] matrix
% Impulse response
xdc.fc = 8e6;  % Center frequency
xdc.bw = .8;  % Fractional bandwidth
tc = gauspuls('cutoff', xdc.fc, xdc.bw, -6, -40);
ti = -tc:1/fs:tc;
xdc.impulseResponse = gauspuls(ti, xdc.fc, xdc.bw);
% Transmit sequence
% We will use a transmit sequence of focused waves, spaced to achieve Nyquist sampling at the elevation focus. The resolution is defined as 
% such that the Nyquist rate should be . For demonstration purposes, we will intentionally undersample the beam spacing in azimuth at .
seq.dx = wl * mean(limz) / xdc.Dx;            % Beam spacing in x
seq.limx = limx;                              % Scan across the whole FOV
seq.x = seq.limx(1):seq.dx:seq.limx(2);       % Place lines in FOV
seq.x = seq.x - mean(seq.x) + mean(seq.limx); % Center the lines
seq.y = 0 * seq.x;                            % Always at y=0
seq.z = 0 * seq.x + mean(limz);               % Sequence focused at middle depth of FOV
seq.foc = [seq.x; seq.y; seq.z];              % Beam focus coordinates
seq.ori = [seq.x; seq.y; seq.z * 0];          % Beam origins are directly above foci (linear scan)
seq.n = length(seq.x);                        % Total number of transmits
% Compute the transmit delays and apodizations for this transmit sequence.
txdel = zeros(xdc.n, seq.n);
txapo = zeros(xdc.n, seq.n);
time0 = zeros(1, seq.n);
for i = 1:seq.n
    % Compute time-of-flight from each element to transmit focus
    txdel(:,i) = sqrt(sum((seq.foc(:,i) - xdc.pos).^2, 1)) / c0;
    txdel(:,i) = txdel(:,i) - min(txdel(:,i));
    % Compute transmit apodization based on beam origin and tx f-number
    ctrs_i = seq.foc(:,i) - xdc.pos;  % Relative element positions
    txapo(:,i) = abs(seq.foc(3,i) ./ (2*ctrs_i(1,:))) >= txfn;
    % Compute the "time zero" of the beam simulation
    time0(i) = interp1(xdc.xc,txdel(:,i),seq.x(i),'pchip','extrap');
end
% Also select the excitation pulse to be applied to the transmit aperture.
ncycles = 0;  % Number of sinusoidal cycles
txexc = 1;    % Impulse by default
if ncycles > 0
    txexc = sin(2*pi*(0:f0/fs:ncycles));
end
% Imaging target
% Define the scatterers that we will image. Here, we choose random point targets throughout the field of view.
nscat = 15;
randvec = @(n, lims) rand(n,1) * (max(lims) - min(lims)) + min(lims);  % lambda function to make random vector
buffer = [+1e-3, -1e-3];  % Buffer to keep points inside the field of view
pos = cat(2, randvec(nscat, limx+buffer), randvec(nscat, limy), randvec(nscat, limz+buffer));
amp = ones(nscat, 1);
% Execute Field II simulation
% Make sure Field II is on your path before starting!
field_init(-1);  % -1 suppresses output
set_sampling(fs);  % Set sampling frequency
set_field('show_times', 0);  % Suppress timing output
% The peak of the pulse-echo waveform has a time offset due to the transmit impulse response, receive impulse response, and the excitation pulse. Determine the time shift according to the centroid of the envelope of the pulse-echo waveform.
% Define the time shift based on the centroid of the total waveform.
wvfm = conv(conv(xdc.impulseResponse(:), xdc.impulseResponse(:)), txexc(:));
wvfm_env = abs(hilbert(wvfm));  % Envelope of waveform
tshift = mean((1:length(wvfm_env))' .* wvfm_env) / mean(wvfm_env) / fs;
% Define transmit and receive apertures based on the presented configuration.
nsubx = 1; nsuby = 5;  % Number of mathematical elements per element
initfoc = [0 0 1];  % Initial electronic focus (will be overwritten)

if xdc.ny == 1 && isinf(xdc.radius)  % If using a 1D linear array
    if isinf(xdc.elevfoc)  % If no elevation lens
        Tx = xdc_linear_array(xdc.nx,xdc.width_x,xdc.width_y,xdc.kerf_x,nsubx,nsuby,initfoc);
        Rx = xdc_linear_array(xdc.nx,xdc.width_x,xdc.width_y,xdc.kerf_x,nsubx,nsuby,initfoc);    
    else  % If elevation lens, use xdc_focused_array
        Tx = xdc_focused_array(xdc.nx,xdc.width_x,xdc.width_y,xdc.kerf_x,xdc.elevfoc,nsubx,nsuby,initfoc);
        Rx = xdc_focused_array(xdc.nx,xdc.width_x,xdc.width_y,xdc.kerf_x,xdc.elevfoc,nsubx,nsuby,initfoc);
        % Account for mechanical lens offset in tshift
        tshift = tshift + (sqrt((xdc.width_y/2)^2 + xdc.elevfoc^2) - xdc.elevfoc) / c0;
        tshift = tshift + (sqrt((xdc.width_y/2)^2 + xdc.elevfoc^2) - xdc.elevfoc) / c0;
    end
elseif xdc.ny == 1 && xdc.radius < inf  % If using a 1D curvilinear array
    if isinf(xdc.elevfoc)  % If no elevation lens
        Tx = xdc_convex_array(xdc.nx,xdc.width_x,xdc.width_y,xdc.kerf_x,xdc.radius,nsubx,nsuby,initfoc);
        Rx = xdc_convex_array(xdc.nx,xdc.width_x,xdc.width_y,xdc.kerf_x,xdc.radius,nsubx,nsuby,initfoc);    
    else  % If elevation lens, use xdc_focused_array
        Tx = xdc_convex_focused_array(xdc.nx,xdc.width_x,xdc.width_y,xdc.kerf_x,xdc.radius,xdc.elevfoc,nsubx,nsuby,initfoc);
        Rx = xdc_convex_focused_array(xdc.nx,xdc.width_x,xdc.width_y,xdc.kerf_x,xdc.radius,xdc.elevfoc,nsubx,nsuby,initfoc);
        % Account for mechanical lens offset in tshift
        tshift = tshift + (sqrt((xdc.width_y/2)^2 + xdc.elevfoc^2) - xdc.elevfoc) / c0;
        tshift = tshift + (sqrt((xdc.width_y/2)^2 + xdc.elevfoc^2) - xdc.elevfoc) / c0;
    end
else  % If using a 2D array, ignore elevation lens and use xdc_2d_array
    enabled = ones(xdc.nx, xdc.ny);
    Tx = xdc_2d_array(xdc.nx,xdc.ny,xdc.width_x,xdc.width_y,xdc.kerf_x,xdc.kerf_y,enabled,nsubx,nsuby,initfoc);
    Rx = xdc_2d_array(xdc.nx,xdc.ny,xdc.width_x,xdc.width_y,xdc.kerf_x,xdc.kerf_y,enabled,nsubx,nsuby,initfoc);
end
% Execute the Field II simulation! We use calc_scat_multi to obtain the channel data.
rfs = cell(seq.n, 1);
samp1s = zeros(seq.n, 1);  % The first sample for each beam in global time
sampNs = zeros(seq.n, 1);  % The last sample for each beam in global time
for i = 1:seq.n
    xdc_impulse(Tx, xdc.impulseResponse);
    xdc_impulse(Rx, xdc.impulseResponse);
    xdc_excitation(Tx, txexc);
    xdc_times_focus(Tx, -1 ,-txdel(:,i)' + time0(i));
    xdc_times_focus(Rx, -1, -tshift * ones(1,xdc.n));
    xdc_apodization(Tx, -1, txapo(:,i)');
    [rfs{i}, t0] = calc_scat_multi(Tx, Rx, pos, amp);
    % Field II seems to work in integer samples
    samp1s(i) = round(t0 * fs);
    sampNs(i) = round(samp1s(i) + size(rfs{i},1) - 1);
end
% Each simulation has its own start time and end time. Align everything to the same time vector.
% Convert to samples
samp1 = min(samp1s);  % Overall first sample
sampN = max(sampNs);  % Overall last sample
nsamps = sampN - samp1 + 1;  % Total number of samples
% Perform alignment
rfdata = zeros(nsamps, xdc.n, seq.n);
for i = 1:seq.n
    rfdata((samp1s(i):sampNs(i)) - samp1 + 1,:,i) = rfs{i};
end
clear rfs
% Set the new aligned time zero to the focal point
time0 = 0 * time0 - samp1/fs + sqrt(sum((seq.foc - seq.ori).^2, 1)) / c0;
% Field II simulations have extremely small floating point values.
% Multiply by a fixed factor to avoid numerical precision errors.
rfdata = rfdata * 1e24;
% Convert to single precision
rfdata = single(rfdata);
% Permute so that dimensions are [samples, transmits, channels]
rfdata = permute(rfdata, [1 3 2]);
% Demodulate and create iqdata
iqdata = rfdata;
for i = 1:seq.n, iqdata(:,i,:) = hilbert(rfdata(:,i,:)); end
iqdata = iqdata .* exp(2j*pi*f0*(0:nsamps-1)'/fs);

% Perform Image Reconstruction (Traditional)
% First, we perform traditional line-by-line reconstruction. In rtbf, this corresponds to the BLKDIAG transmit beamforming mode. It is called "block diagonal" because the transmit beamformer is effectively a block-diagonal matrix (see D Hyun, MICCAI 2022 ASMUS workshop).
% Select pixel grid
% Define a desired pixel grid. Here, we will create pixels coming straight down from the beam origins.
xi = seq.x;
yi = 0;
zi = limz(1):wl/2:limz(2);  % wl/2 pixel spacing
if mod(length(zi),2) == 0, zi = zi(1:end-1); end
zi = zi - mean(zi) + mean(limz);
[zz, xx, yy] = ndgrid(zi, xi, yi);
% Pixel positions should have (x,y,z) in first dimension
pxpos = cat(4, xx, yy, zz); pxpos = permute(pxpos, [4 1 2 3]);
Prepare MATLAB structs to input into rtbf
Define the Focus struct.
F.operator = 'Focus';
F.pxpos = single(pxpos);  % Should have size [3, ...]
F.elpos = single(xdc.pos);  % Should have size [3, nelems]
F.txdat = single(seq.foc);  % Should have size [3, nxmits]
F.fs = fs; F.fc = 0; F.fdemod = f0; F.c0 = c0;
F.txbf_mode = 'BLKDIAG';  % Traditional line-by-line beamforming
F.rxbf_mode = 'IDENTITY';  % Do not sum the receives; leave them as is.
F.rxfnum = rxfn;
F.time0 = single(time0);
% (Instead of plane wave focusing, you can also perform virtual source focusing by supplying F.txdat = virtual_source_positions; as a [3, nxmits] array.)
% Define the ChannelSum struct.
C.operator = 'ChannelSum';
C.axis = -1;  % Sum the last dimension. (Default, so this line could be omitted.)
% Define the Hilbert struct. (No extra parameters.)
H.operator = 'Hilbert';
% Define the Bmode struct. Options include logarithmic and power compression.
B.operator = 'Bmode';
B.compression = 'log';  % Alternatively 'power', which has parameter B.gamma.
% Reconstruct
% Perform traditional transmit beamforming.
clear rtbf_mex
bimg1 = rtbf_mex(iqdata, F, C, B);
imagesc(xi*1e3, zi*1e3, bimg1, [-60 0]+max(bimg1(:)));
axis image; colorbar
title('Traditional Tx Beamforming')

% Perform Image Reconstruction (Same grid, retrospective transmit beamforming)
% This time, let's perform reconstruction with retrospective transmit beamforming (RTB) for a dynamic transmit focus. In rtbf, this corresponds to the 'SUM' transmit beamforming mode, where all contributing transmit events are included according to the specified transmit f-number.
% Use Focus for transmit aperture synthesis
% The Focus struct is mostly the same, but with two key differences.
F2 = F;
F2.txbf_mode = 'SUM';  % This time, sum across transmits as well
F2.txfnum = txfn;  % This time, specify virtual source f-number
% Reconstruct
% Perform RTB. This method treats the transmit foci as virtual sources and beamforms accordingly.
clear rtbf_mex
bimg2 = rtbf_mex(iqdata, F2, C, B);
imagesc(xi*1e3, zi*1e3, bimg2, [-60 0]+max(bimg2(:)));
axis image; colorbar; grid on
title('Retrospective Tx Beamforming')

% Perform Image Reconstruction (Arbitrary grid)
% With RTB, our reconstruction is no longer bound by the transmit "beams". Transmit and receive beamforming can be performed at arbitrary pixel positions.
% Select a finer pixel grid
% Let's select a pixel grid with finer isotropic sampling of /3 pixel spacing.
dx = wl/3; dz = wl/3;
xi = limx(1):dx:limx(2);
xi = xi - mean(xi) + mean(limx);
yi = 0;
zi = limz(1):dz:limz(2);  % wl/2 pixel spacing
zi = zi - mean(zi) + mean(limz);
[zz, xx, yy] = ndgrid(zi, xi, yi);
% Pixel positions should have (x,y,z) in first dimension
pxpos = cat(4, xx, yy, zz); pxpos = permute(pxpos, [4 1 2 3]);
% Use Focus to specify new pixel grid
% We can re-use the Focus struct above, replacing only the pixel positions.
F3 = F2;
F3.pxpos = single(pxpos);  % This time, change the pixel positions
% Reconstruct
% Perform RTB on this finer pixel grid. Observe the improved visual quality of the image. (Note: once dx and dz are reduced beyond the Nyquist rate, there is no benefit to finer sampling. The same results can be achieved by simply interpolating the final image.)
clear rtbf_mex
bimg3 = rtbf_mex(iqdata, F3, C, B);
imagesc(xi*1e3, zi*1e3, bimg3, [-60 0]+max(bimg3(:)));
axis image; colorbar; grid on
title('Retro. Tx. BF with a Finer Grid')

% Perform Image Reconstruction (Refocus)
% Finally, let us perform RTB using the "refocus" technique (N Bottenus, TUFFC 2018 and R Ali et al., TUFFC 2019). Refocus does not rely on the virtual source model, and instead relies on pseudoinverting the transmit sequence directly to recover the multistatic (i.e. single-element transmit) synthetic aperture dataset (D Hyun, MICCAI 2022). This code can be used to beamform arbitrary transmit sequences in real-time (D Hyun et al., IUS 2021). The only inputs needed are the transmit delays, apodizations, and sampling frequency.
% Create a Refocus struct
% Define the Refocus struct.
R.operator = 'Refocus';
R.txdel = single(-txdel);  % Opposite sign convention with Field II
R.txapo = single(txapo);
R.fs = fs;
% Modify Focus to work with refocused data
% The Focus struct is mostly the same, but with two key differences.
F4.operator = 'Focus';
F4.pxpos = single(pxpos);
F4.elpos = single(xdc.pos);
F4.txdat = single(xdc.pos);  % This time, change the pixel positions
F4.fs = fs; F4.fc = f0; F4.c0 = c0;
F4.time0 = zeros(xdc.n, 1, 'single') - samp1 / fs;
F4.txbf_mode = 'SUM';  % Sum the transmits.
F4.rxbf_mode = 'IDENTITY';  % Do not sum the receives; leave them as is.
F4.txfnum = txfn; F4.rxfnum = rxfn;
% Reconstruct
% Perform RTB with refocus. Observe that the point targets away from the elevation focus (20 mm) have improved focusing.
clear rtbf_mex
bimg4 = rtbf_mex(rfdata, H, R, F4, C, B);
imagesc(xi*1e3, zi*1e3, bimg4, [-60 0]+max(bimg4(:)));
axis image; colorbar; grid on
title('Refocused Tx Beamforming')
