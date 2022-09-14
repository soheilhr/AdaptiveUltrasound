%% run all the files
clear all; close all

myFolder = 'C:/Users/Shrey/Documents/Stanford/Ultrasound/EchoNetFrames/0X1A0A263B22CCD966/frames';
% Check if the folder exists
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: the following folder does not exist:\n%s\n Please specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new folder
    if myFolder == 0 % user clicked cancel
        return;
    end
end

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.tif'); % specifies which type of files we want from the folder in case there is something other than .tif present there
theFiles = dir(filePattern); % directory: all the frames for this patient
for k = 1 : 2
    % Get and read the file
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    frame_mat = imread(fullFileName);
    frame_mat_bw = rgb2gray(frame_mat);
    frame_mat_bw_mask = imbinarize(frame_mat_bw, 0.5); % threshold: set to 0.5, can consider otsu threshold or other

    % Save the black and white version of the frame and the mask of the
    % frame
    save('framebw.mat', 'frame_mat_bw');
    save('framebwmask.mat', 'frame_mat_bw_mask');

    % Define the x and z axes for the template and phantom image
    x_axis = linspace(-.026, .026, 112);
    z_axis = linspace(1.5e-4, .0444299897172237, 112);

    % Define the Bmode, sector, x axis, and z axis in the form of a struct
    img.bmode = frame_mat_bw;
    img.sector = frame_mat_bw_mask;
    img.x_axis = x_axis;
    img.z_axis = z_axis;

    %-- define some useful parameters
    bckDensity = 20;

    %-- display B mode image
    figure; imagesc(x_axis*1e3,z_axis*1e3,img.bmode); 
    axis image; colormap(gray); colorbar;
    % hold on; contour(x_axis*1e3,z_axis*1e3,img.sector,[0.5 0.5],'y','linewidth',3);
    xlabel('mm'); ylabel('mm'); title('Template');
    bmode_folder = 'C:\Users\Shrey\Documents\Stanford\Ultrasound\EchoNetFrames\0X1A0A263B22CCD966\bmode';
    saveas(gcf, fullfile(bmode_folder, baseFileName), 'jpg');

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
    fc = 2e6;
    lambda = 1540/fc;
    psf_length = 3 * lambda / 2; % point spread function length
    psf_width = 1.206 * lambda * f_number;        % point spread function width
    N_sca_speckle = round( bckDensity * height * width / psf_length / psf_width );            
    
    %-- x positions of the scatterers
    xx = random('unif',x_min,x_max,N_sca_speckle,1);

    %-- z positions of the scatterers
    zz = random('unif',z_min,z_max,N_sca_speckle,1);  
    [xm,zm] = meshgrid(x_axis,z_axis);
    %--amplitudes of the scatterers computed from the intensity of the template image
    ampl = interp2(xm,zm,double(img.bmode),xx,zz,'linear');
    I = ampl;
    dyn_range = 10; % can change dynamic range 
    exp = (dyn_range/20)*(I/max(ampl) - 1);
    ampl = 10.^exp; % lines 78 - 79 are the formula given in the MUST lecture
                        
    %-- select only scatterers inside the region of interest
    start_depth = lambda; % 7.7E-4
    end_depth = z_max; % .0444
    IN = hypot(xx,zz)>=lambda & hypot(xx,zz)<=end_depth;   % determines the region of interest that is valid -- where xx+zz is more than lambda and less than end depth
    % hypot: computes the square root of the sum of the squares of 2 inputs
    % (xx, zz)
    IN = IN & abs((angle(xx+1i*zz)-pi/2))<50/180*pi; % angle(z): returns the phase angle in the interval [-pi, pi] for each element of a complex array z        
    sca_x = xx(IN); % scatterers in x
    sca_z = zz(IN); % scatterers in z
    ampl = ampl(IN); % amplitude of the specific point
    
    %-- save results
    %save('../fig_1_temp.mat','sca_x','sca_z','ampl', 'img.bmode','x_axis','z_axis','img.sector');        
    

    %-- display scatterers that will be used for the simulation
    figure; scatter(sca_x*1e3,sca_z*1e3,1,ampl,'linewidth',3);
    colorbar; grid on;
    axis equal; axis ij;
    xlabel('x [mm]');
    ylabel('z [mm]');
    title('scatterers');
    set(gca,'fontsize',16);
    
    phantom_folder = 'C:\Users\Shrey\Documents\Stanford\Ultrasound\EchoNetFrames\0X1A0A263B22CCD966\phantom';
    saveas(gcf, fullfile(phantom_folder, baseFileName), 'jpg');

end