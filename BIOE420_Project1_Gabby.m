

%% BIOE420 Project 1: CT Image Reconstruction 
% Group Members: Lillian Myers, Julie Karasik, Kenneth Missig, Katherine
% Chan, Gabby Zanga

%% A.
% Exercise 1.
clear; clc;

%Extraction of files

files = dir(fullfile('Crystal 360 Images','*.bmp'));

% sort image files by angle 
angles = zeros(length(files), 1);
for k = 1:length(files)

    name = files(k).name;
    token = regexp(name, '(\d+\.?\d*)\.bmp$', 'tokens');
    angles(k) = str2double(token{1}{1});

end

[angles, order] = sort(angles);

files = files(order);

% use only angles 0-359, since 0 degrees and 360 degrees are the same
use = (angles >= 0) & (angles < 360);

files = files(use);
angles = angles(use);

% read first image to determine size
    img = imread(fullfile('Crystal 360 Images', files(1).name));

    % Convert to double
    img = double(img);

    % Get image dimensions
    [Ny, Nx] = size(img);

% load images and create sinogram
    x0 = 500;
    y0 = 750;

    angleStep = 2;
    Nproj = 180;
    sinogram = zeros(Nproj, Nx);

for k = 1:Nproj

    % Load image
    img = imread(fullfile('Crystal 360 Images', files(k).name));

    % Convert to double
    img = double(img);

    % Extract horizontal line at y = 750
    sinogram(k,:) = img(y0,:);

end


% display sonogram
    figure;

    imagesc(angles, 1:Nx, sinogram');

    colormap gray;
    axis image;
    colorbar;
    
    xlabel('Angle (degrees)');
    ylabel('Detector position');
    
    title('Sinogram for Crystal 360 Dataset');


    
% Excise 2. 


% Excise 3.


% Excise 4. 


%% B. 
% Excise 1.
clear; clc;


% Excise 2.


% Excise 3.


% Excise 4.


% Excise 5. 
