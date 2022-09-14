% Extract data from videos from EchoNET
outputFolder = '/home/shreyavn/Ultrasound/EchoNetProcessing/frames_test_patient/EchoNetFrames/0X1A2A76BDB5B98BED/frames/';
obj = VideoReader('0X1A2A76BDB5B98BED.avi');
vid = read(obj);

frames = obj.NumFrames;

for x = 1:frames
    imwrite(vid(:,:,:,x),strcat(outputFolder,'frame-', num2str(x), '.tif'));
%     outputBaseFileName = sprintf('Frame %4.4', frame);
%     outputFileName = fullfile(outputFolder, outputBaseFileName);
%     imwrite(mat_file, outputFileName, 'mat')
end

% TODO: make it so that we can run for all files and save each set of
% associated frames to its own separate folder. 

% %% for now, convert one of the frames to a mat file: 
% 
% I = zeros(112,112,3);
% frame_1_mat = imread('frame-1.tif');
% frame_1_mat_bw = rgb2gray(frame_1_mat);
% % threshold: set to 0.5, can consider otsu threshold or other
% frame_1_mat_bw_mask = imbinarize(frame_1_mat_bw, 0.5);
% save('frame1bw.mat', 'frame_1_mat_bw');
% save('frame1bwmask.mat', 'frame_1_mat_bw_mask');