%% BIOE420 Project 1: CT Image Reconstruction 
% Group Members: Lillian Myers, Julie Karasik, Kenneth Missig, Katherine
% Chan


%% BIOE420 Project 1
% Exercises 1 and 2
% Manual Filtered Back Projection
clear; clc; close all;

%% =========================================================
% EXERCISE 1: CRYSTAL DATASET
% ==========================================================

% Unzip crystal dataset
unzip('Data_Crystal_360.zip', 'Crystal 360 Images');

% Find all BMP files, including subfolders
files = dir(fullfile('Crystal 360 Images', '**', '*.bmp'));

if isempty(files)
    error('No BMP images found for Crystal dataset.');
end

%% Sort images by projection angle

angles = zeros(length(files), 1);

for k = 1:length(files)

    token = regexp(files(k).name, ...
        '(\d+\.?\d*)\.bmp$', 'tokens');

    if isempty(token)
        error('Cannot extract angle from %s', files(k).name);
    end

    angles(k) = str2double(token{1}{1});

end

[angles, order] = sort(angles);
files = files(order);

% Keep all angles from 0 to 358 degrees
use = (angles >= 0) & (angles < 360);

angles = angles(use);
files = files(use);

%% Read first image

img = imread(fullfile(files(1).folder, files(1).name));

if size(img, 3) == 3
    img = rgb2gray(img);
end

img = double(img);

[Ny, Nx] = size(img);

% Center position specified in assignment
x0 = 500;
y0 = 750;

if x0 > Nx || y0 > Ny
    error('Crystal center position exceeds image dimensions.');
end

%% Create Crystal Sinogram (0-358 degrees)

Nproj = length(angles);
sinogram = zeros(Nx, Nproj);

for k = 1:Nproj

    img = imread(fullfile(files(k).folder, files(k).name));

    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    img = double(img);

    % Extract horizontal line at y = 750
    sinogram(:, k) = img(y0, :)';

end

%% Display Full Crystal Sinogram (0-360 degrees)

figure;

imagesc(angles, 1:Nx, sinogram);

colormap gray;
axis normal;
colorbar;

xlim([0 360]);
xticks(0:45:360);

xlabel('Projection Angle (degrees)');
ylabel('Detector Position (pixels)');
title('Crystal Sinogram (0-360 Degrees)');

%% Select 0-178 Degrees for Filtered Back Projection

% Parallel-beam reconstruction needs only 180 degrees
use180 = (angles >= 0) & (angles < 180);

theta = angles(use180);
projectionData = sinogram(:, use180);

NprojFBP = length(theta);

%% Apply Ram-Lak Filter - Crystal

Nfft = 2^nextpow2(2 * Nx);

freq = (0:Nfft-1)';
freq(freq > Nfft/2) = freq(freq > Nfft/2) - Nfft;

rampFilter = abs(freq) / Nfft;

filteredSinogram = zeros(Nx, NprojFBP);

for k = 1:NprojFBP

    projectionFFT = fft(projectionData(:, k), Nfft);

    filteredFFT = projectionFFT .* rampFilter;

    filteredProjection = real(ifft(filteredFFT));

    filteredSinogram(:, k) = filteredProjection(1:Nx);

end

%% Manual Back Projection - Crystal

N = Nx;

[x, y] = meshgrid((1:N)-x0, (1:N)-x0);

reconstructedCrystal = zeros(N, N);

detectorPositions = (1:Nx) - x0;

for k = 1:NprojFBP

    angle = deg2rad(theta(k));

    t = x*cos(angle) + y*sin(angle);

    backProjection = interp1( ...
        detectorPositions, ...
        filteredSinogram(:, k), ...
        t, 'linear', 0);

    reconstructedCrystal = reconstructedCrystal + backProjection;

end

if NprojFBP > 1
    deltaTheta = deg2rad(mean(diff(theta)));
    reconstructedCrystal = reconstructedCrystal * deltaTheta;
end

%% Display Crystal Reconstruction

figure;

imagesc(reconstructedCrystal);

colormap gray;
axis image;
colorbar;

xlabel('X Position (pixels)');
ylabel('Y Position (pixels)');
title('Crystal Cross-Section - Filtered Back Projection');

%% Display Crystal Results Side by Side

figure;

subplot(1,2,1);

imagesc(angles, 1:Nx, sinogram);
colormap gray;
axis normal;
colorbar;

xlim([0 360]);
xticks(0:45:360);

xlabel('Angle (degrees)');
ylabel('Detector Position');
title('Crystal Sinogram (0-360 Degrees)');

subplot(1,2,2);

imagesc(reconstructedCrystal);
colormap gray;
axis image;
colorbar;

xlabel('X Position');
ylabel('Y Position');
title('Crystal Reconstructed Slice');


%% Step 1: Load Insect Head Dataset

zipFile = 'Dataset_Insect_Head_360.zip';
outputFolder = 'Insect Head Images';

% Only unzip if the folder does not already exist
if ~isfolder(outputFolder)
    unzip(zipFile, outputFolder);
end

% Find TIFF images in all subfolders
files = [ ...
    dir(fullfile(outputFolder, '**', '*.tif')); ...
    dir(fullfile(outputFolder, '**', '*.tiff'))];

if isempty(files)
    error('No TIFF images found. Check the dataset folder.');
end

fprintf('Found %d TIFF files.\n', length(files));

%% Step 2: Extract Angles from Filenames

angles = [];
validFiles = files([]);

for k = 1:length(files)

    name = files(k).name;

    % Handles filenames such as:
    % 202009231201PM49_ -0.0.tif
    % 202009231202PM05_ 2.0 - Copy.tif
    % 202009231202PM05_ 4.0.tif

    token = regexp(name, ...
        '_\s*([+-]?\d+\.?\d*)\s*(?:-\s*Copy)?\.tiff?$', ...
        'tokens', 'once', 'ignorecase');

    if isempty(token)
        fprintf('Skipping invalid filename: %s\n', name);
        continue;
    end

    angle = str2double(token{1});

    if ~isfinite(angle)
        fprintf('Skipping invalid angle: %s\n', name);
        continue;
    end

    % Convert negative angles to 0-360 degrees
    angle = mod(angle, 360);

    angles(end+1, 1) = angle;
    validFiles(end+1, 1) = files(k);

end

if isempty(validFiles)
    error('No valid projection filenames found.');
end

%% Step 3: Sort Images and Remove Duplicate Angles

[angles, order] = sort(angles);
files = validFiles(order);

% Prefer original files over "- Copy" files
isCopy = contains({files.name}, '- Copy', 'IgnoreCase', true);
[~, priorityOrder] = sort(isCopy);

angles = angles(priorityOrder);
files = files(priorityOrder);

% Remove duplicate angles
[angles, uniqueIdx] = unique(angles, 'stable');
files = files(uniqueIdx);

% Sort again by angle
[angles, order] = sort(angles);
files = files(order);

fprintf('Found %d unique projection angles.\n', length(angles));
fprintf('Angle range: %.1f to %.1f degrees.\n', ...
    min(angles), max(angles));

%% Step 4: Find First Readable Image

firstValid = [];

for k = 1:length(files)

    try
        img = imread(fullfile(files(k).folder, files(k).name));

        if size(img, 3) == 3
            img = rgb2gray(img);
        end

        firstValid = k;
        break;

    catch ME
        fprintf('Cannot read: %s\n', files(k).name);
        fprintf('Reason: %s\n', ME.message);
    end

end

if isempty(firstValid)
    error('No readable TIFF images found.');
end

img = double(img);

[Ny, Nx] = size(img);

fprintf('Image dimensions: %d x %d pixels.\n', Ny, Nx);

%% Step 5: Define Cross-Section Position

% Use center of insect head image
x0 = round(Nx/2);
y0 = round(Ny/2);

fprintf('Selected center: (%d, %d)\n', x0, y0);

%% Step 6: Create Sinogram

Nproj = length(angles);

% Rows = detector positions
% Columns = projection angles
sinogram = zeros(Nx, Nproj);

validProjection = false(1, Nproj);

for k = 1:Nproj

    filename = fullfile(files(k).folder, files(k).name);

    try

        img = imread(filename);

        if size(img, 3) == 3
            img = rgb2gray(img);
        end

        img = double(img);

        % Verify dimensions
        if size(img,1) ~= Ny || size(img,2) ~= Nx
            error('Image dimensions do not match.');
        end

        % Extract horizontal line through center
        sinogram(:,k) = img(y0,:)';

        validProjection(k) = true;

    catch ME

        fprintf('Skipping unreadable image: %s\n', ...
            files(k).name);
        fprintf('Reason: %s\n', ME.message);

    end

end

% Remove unreadable projections
sinogram = sinogram(:,validProjection);
angles = angles(validProjection);

Nproj = length(angles);

if Nproj < 2
    error('Not enough readable projections.');
end

fprintf('Successfully loaded %d projections.\n', Nproj);

%% Step 7: Display Full Sinogram (0-360 Degrees)

figure;

imagesc(angles, 1:Nx, sinogram);

colormap gray;
axis normal;
colorbar;

xlim([0 360]);
xticks(0:45:360);

xlabel('Projection Angle (degrees)');
ylabel('Detector Position (pixels)');
title('Insect Head Sinogram (0-360 Degrees)');

%% Step 8: Select Angles for Reconstruction

% Use 0-180 degrees for parallel-beam FBP
use180 = (angles >= 0) & (angles < 180);

theta = angles(use180);
projectionData = sinogram(:,use180);

NprojFBP = length(theta);

if NprojFBP < 2
    error('Not enough projections between 0 and 180 degrees.');
end

%% Step 9: Apply Ram-Lak Filter Manually

% Zero-padding for FFT filtering
Nfft = 2^nextpow2(2*Nx);

% Frequency indices in FFT order
freq = (0:Nfft-1)';
freq(freq > Nfft/2) = freq(freq > Nfft/2)-Nfft;

% Ram-Lak frequency response
rampFilter = abs(freq)/Nfft;

% Filtered projections
filteredSinogram = zeros(Nx,NprojFBP);

for k = 1:NprojFBP

    % Fourier transform
    projectionFFT = fft(projectionData(:,k),Nfft);

    % Apply Ram-Lak filter
    filteredFFT = projectionFFT .* rampFilter;

    % Inverse Fourier transform
    filteredProjection = real(ifft(filteredFFT));

    % Keep original detector length
    filteredSinogram(:,k) = filteredProjection(1:Nx);

end

%% Step 10: Manual Back Projection

% Reconstruction grid
N = Nx;

[x,y] = meshgrid((1:N)-x0,(1:N)-x0);

% Initialize reconstructed image
reconstructedInsect = zeros(N,N);

% Detector positions relative to rotation center
detectorPositions = (1:Nx)-x0;

% Angular integration weights
edges = [0; ...
    (theta(1:end-1)+theta(2:end))/2; ...
    180];

weights = deg2rad(diff(edges));

for k = 1:NprojFBP

    angle = deg2rad(theta(k));

    % Detector coordinate for every pixel
    t = x*cos(angle)+y*sin(angle);

    % Interpolate filtered projection
    backProjection = interp1( ...
        detectorPositions, ...
        filteredSinogram(:,k), ...
        t, 'linear', 0);

    % Add weighted contribution
    reconstructedInsect = reconstructedInsect + ...
        backProjection*weights(k);

end

%% Step 11: Display Reconstructed Cross-Section

figure;

imagesc(reconstructedInsect);

colormap gray;
axis image;
colorbar;

xlabel('X Position (pixels)');
ylabel('Y Position (pixels)');
title('Insect Head - Manual Filtered Back Projection');

%% Step 12: Display Both Results Side by Side

figure;

subplot(1,2,1);

imagesc(angles,1:Nx,sinogram);

colormap gray;
axis normal;
colorbar;

xlim([0 360]);
xticks(0:45:360);

xlabel('Projection Angle (degrees)');
ylabel('Detector Position (pixels)');
title('Insect Head Sinogram (0-360 Degrees)');

subplot(1,2,2);

imagesc(reconstructedInsect);

colormap gray;
axis image;
colorbar;

xlabel('X Position (pixels)');
ylabel('Y Position (pixels)');
title('Insect Head Reconstructed CT Slice');







