%% for now, convert one of the frames to a mat file: 
clear all; close all

frame_1_mat = imread('frame-1.tif');
frame_1_mat_bw = rgb2gray(frame_1_mat);
% threshold: set to 0.5, can consider otsu threshold or other
frame_1_mat_bw_mask = imbinarize(frame_1_mat_bw, 0.5);
% save('frame1bw.mat', 'frame_1_mat_bw');
save('frame1bwmask.mat', 'frame_1_mat_bw_mask');


% x_axis = [-.026:.000465:.026];
x_axis  = linspace(-.075, .075, 112);
% x_axis = linspace(-.039, .033, 112);
% x_axis = linspace(-.030, .035, 112);
x_axis_template = linspace(-.075, .075, 112);
z_axis = linspace(0, .15, 112);
%z_axis = linspace(1.5e-4, .0444299897172237, 112);
% z_axis = [1.50e-4:.000396:0.0444299897172237];

img.bmode = frame_1_mat_bw;
img.sector = frame_1_mat_bw_mask;
img.x_axis = x_axis;
img.z_axis = z_axis;

%-- define some useful parameters
bckDensity = 5;

%-- load template data
%load('../data/template.mat'); % 1 x 1 struct with 4 fields (bmode, sector mask, x axis, z axis)
%-- display template image
figure; imagesc(x_axis*1e3,z_axis*1e3,img.bmode); 
axis image; colormap(gray); colorbar;
% hold on; contour(x_axis*1e3,z_axis*1e3,img.sector,[0.5 0.5],'y','linewidth',3);
xlabel('mm'); ylabel('mm'); title('Template');

%-- define a medium from a set of point scatterers ramdomly distributed    
x_min = min(x_axis);
x_max = max(x_axis);
z_min = min(z_axis);
z_max = max(z_axis);
width = x_max-x_min;
height = z_max-z_min;

%-- select the number of scatterers 
%-- the goal is to have "bckDensity" of scatterers per resolution cell
f_number = 1.;
fc = 5e6;
lambda = 1540/fc;
psf_length = 3 * lambda / 2; % point spread function length; originally 3
psf_width = 1.206 * lambda * f_number;        % point spread function width; orig 1.206
N_sca_speckle = round( bckDensity * height * width / psf_length / psf_width );            
%-- x positions of the scatterers
xx = random('unif',x_min,x_max,N_sca_speckle,1);
%-- z positions of the scatterers
zz = random('unif',z_min,z_max,N_sca_speckle,1);        
[xm,zm] = meshgrid(x_axis,z_axis);
%--amplitudes of the scatterers computed from the intensity of the template image
ampl = interp2(xm,zm,double(img.bmode),xx,zz,'linear');
I = ampl;

%-- Insert your own code based on equation (1) from the hands-on description
% dyn_range = 10;
% exp = (dyn_range/20)*(I/max(ampl) - 1);
% ampl = 10.^exp;
%                 
%-- select only scatterers inside the region of interest
start_depth = lambda;
end_depth = z_max;
IN =  zz > .00088524 & zz >= -1.0088*xx + .01635 & zz>= 1.04497*xx + .00135;
% IN = hypot(xx,zz)>=lambda & hypot(xx,zz)<=end_depth;  
% math it out -- dont use line above and below
% subtract off something from z and angle with x; derive myself
% IN = IN & abs((angle(xx+1i*zz)-pi/2));%<60/180*pi;        
sca_x = xx(IN);
sca_z = zz(IN);        
ampl = ampl(IN);

%-- save results
save('../test_phantom_patient1_bd5_new_axis.mat','sca_x','sca_z','ampl','x_axis','z_axis');        

%-- display scatterers that will be used for the simulation
figure; scatter(sca_x*1e3,sca_z*1e3,1,ampl,'linewidth',3);
colorbar; grid on;
axis equal; axis ij;
xlabel('x [mm]');
ylabel('z [mm]');
title('scatterers');
set(gca,'fontsize',16);
