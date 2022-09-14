% Example of use of the new Field II program running under Matlab
%
% This example shows how flow can simulated
%
% This script assumes that the field_init procedure has been called
%
% Example by Joergen Arendt Jensen, March 22, 2011.
% Generate the transducer apertures for send and receive
f0=3e6; % Transducer center frequency [Hz]
fs=100e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda=c/f0; % Wavelength
element_height=5/1000; % Height of element [m]
kerf=0.1/1000; % Kerf [m]
focus=[0 0 70]/1000; % Fixed focal point [m]
% Generate aperture
aperture = xdc_linear_array (128, lambda/2, element_height, kerf, 1, 1,focus);
% Set the impulse response and excitation of the emit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (aperture, impulse_response);
excitation=sin(2*pi*f0*(0:1/fs:8/f0));
xdc_excitation (aperture, excitation);
% Set the seed of the random number generator
rng('shuffle')
% Initialize the ranges for the scatterers
% Notice that the coordinates are in meters
x_range=0.015; % x range for the scatterers [m]
y_range=0.015; % y range for the scatterers [m]
z_range=0.015; % z range for the scatterers [m]
z_offset=0.70; % Offset of the mid-point of the scatterers [m]
R=0.005; % Radius of blood vessel [m]
% Set the number of scatterers. It should be roughly
% 10 scatterers per resolution cell
c=1540; % Ultrasound propagation velocity [m/s]
f0=3e6; % Center frequency of transducer [Hz]
lambda=c/f0;
N=round(10*x_range/(5*lambda)*y_range/(5*lambda)*z_range/(lambda*2));
disp([num2str(N),'Scatterers'])
% Generate the coordinates and amplitude
% Coordinates are rectangular within the range.
% The amplitude has a Gaussian distribution.
x=x_range*(rand(1,N)-0.5);
y=y_range*(rand(1,N)-0.5);
z=z_range*(rand(1,N)-0.5);
% Find which scatterers that lie within the blood vessel 
r=(y.^2+z.^2).^0.5;
within_vessel= (r < R)';
% Assign an amplitude and a velocity for each scatterer
v0=0.5; % Largest velocity of scatterers [m/s]
velocity=v0*(1-(r/R).^2).*within_vessel';
blood_to_stationary= 0.1; % Ratio between amplitude of blood to stationary tissue
amp=randn(N,1).*((1-within_vessel) + within_vessel*blood_to_stationary);
% Calculate a suitable Tprf
theta=45/180*pi;
f_max=2*v0*cos(theta)/c*f0;
fprf=3*f_max;
Tprf=1/fprf; % Time between pulse emissions [s]
Nshoots=128; % Number of shoots
% Find the response by calling field
for i=1:Nshoots
    % Generate the rotated and offset block of sample
    xnew=x*cos(theta)+z*sin(theta);
    znew=z*cos(theta)-x*sin(theta) + z_offset;
    scatterers=[xnew; y; znew;]' ;
    % Calculate the received response
    [v, t1]=calc_scat(aperture, aperture, scatterers, amp);
    % Store the result
    image_data(1:max(size(v)),i)=v;
    times(i) = t1;
    % Propagate the scatterers and alias them
    % to lie within the correct range
    x=x + velocity*Tprf;
    outside_range= (x > x_range/2);
    x=x - x_range*outside_range;
end
% Here the display of the data is inserted
plot(image_data)