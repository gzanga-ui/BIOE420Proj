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

figure('Position', [100, 100, 1200, 500]);

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

figure('Position', [100, 100, 1400, 600]);

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


%% =========================================================
% EXERCISE 2: INSECT HEAD DATASET
% ==========================================================

% Unzip insect head dataset
unzip('Dataset_Insect_Head_360.zip', 'Insect Head Images');

% Find BMP files in subfolders
files2 = dir(fullfile('Insect Head Images', '**', '*.bmp'));

if isempty(files2)
    error('No BMP images found for Insect Head dataset.');
end

%% Sort images by projection angle

angles2 = zeros(length(files2), 1);

for k = 1:length(files2)

    token = regexp(files2(k).name, ...
        '(\d+\.?\d*)\.bmp$', 'tokens');

    if isempty(token)
        error('Cannot extract angle from %s', files2(k).name);
    end

    angles2(k) = str2double(token{1}{1});

end

[angles2, order2] = sort(angles2);
files2 = files2(order2);

% Keep all angles from 0 to 358 degrees
use2 = (angles2 >= 0) & (angles2 < 360);

angles2 = angles2(use2);
files2 = files2(use2);

%% Read first image

img2 = imread(fullfile(files2(1).folder, files2(1).name));

if size(img2, 3) == 3
    img2 = rgb2gray(img2);
end

img2 = double(img2);

[Ny2, Nx2] = size(img2);

% Center position for insect head
x02 = round(Nx2/2);
y02 = round(Ny2/2);

%% Create Insect Head Sinogram (0-358 degrees)

Nproj2 = length(angles2);
sinogram2 = zeros(Nx2, Nproj2);

for k = 1:Nproj2

    img2 = imread(fullfile(files2(k).folder, files2(k).name));

    if size(img2, 3) == 3
        img2 = rgb2gray(img2);
    end

    img2 = double(img2);

    % Extract horizontal line
    sinogram2(:, k) = img2(y02, :)';

end

%% Display Full Insect Head Sinogram (0-360 degrees)

figure('Position', [100, 100, 1200, 500]);

imagesc(angles2, 1:Nx2, sinogram2);

colormap gray;
axis normal;
colorbar;

xlim([0 360]);
xticks(0:45:360);

xlabel('Projection Angle (degrees)');
ylabel('Detector Position (pixels)');
title('Insect Head Sinogram (0-360 Degrees)');

%% Select 0-178 Degrees for Filtered Back Projection

use180_2 = (angles2 >= 0) & (angles2 < 180);

theta2 = angles2(use180_2);
projectionData2 = sinogram2(:, use180_2);

NprojFBP2 = length(theta2);

%% Apply Ram-Lak Filter - Insect Head

Nfft2 = 2^nextpow2(2 * Nx2);

freq2 = (0:Nfft2-1)';
freq2(freq2 > Nfft2/2) = ...
    freq2(freq2 > Nfft2/2) - Nfft2;

rampFilter2 = abs(freq2) / Nfft2;

filteredSinogram2 = zeros(Nx2, NprojFBP2);

for k = 1:NprojFBP2

    projectionFFT2 = fft(projectionData2(:, k), Nfft2);

    filteredFFT2 = projectionFFT2 .* rampFilter2;

    filteredProjection2 = real(ifft(filteredFFT2));

    filteredSinogram2(:, k) = filteredProjection2(1:Nx2);

end

%% Manual Back Projection - Insect Head

N2 = Nx2;

[x2, y2] = meshgrid((1:N2)-x02, (1:N2)-x02);

reconstructedInsect = zeros(N2, N2);

detectorPositions2 = (1:Nx2) - x02;

for k = 1:NprojFBP2

    angle2 = deg2rad(theta2(k));

    t2 = x2*cos(angle2) + y2*sin(angle2);

    backProjection2 = interp1( ...
        detectorPositions2, ...
        filteredSinogram2(:, k), ...
        t2, 'linear', 0);

    reconstructedInsect = reconstructedInsect + backProjection2;

end

if NprojFBP2 > 1
    deltaTheta2 = deg2rad(mean(diff(theta2)));
    reconstructedInsect = reconstructedInsect * deltaTheta2;
end

%% Display Insect Head Reconstruction

figure;

imagesc(reconstructedInsect);

colormap gray;
axis image;
colorbar;

xlabel('X Position (pixels)');
ylabel('Y Position (pixels)');
title('Insect Head Cross-Section - Filtered Back Projection');

%% Display Insect Head Results Side by Side

figure('Position', [100, 100, 1400, 600]);

subplot(1,2,1);

imagesc(angles2, 1:Nx2, sinogram2);
colormap gray;
axis normal;
colorbar;

xlim([0 360]);
xticks(0:45:360);

xlabel('Angle (degrees)');
ylabel('Detector Position');
title('Insect Head Sinogram (0-360 Degrees)');

subplot(1,2,2);

imagesc(reconstructedInsect);
colormap gray;
axis image;
colorbar;

xlabel('X Position');
ylabel('Y Position');
title('Insect Head Reconstructed Slice');


%% Exercise 3.


%% Exercise 4.












