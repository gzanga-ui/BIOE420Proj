%% BIOE420 Project 1: CT Image Reconstruction
% Group Members: Julie Karasik, Gabrielle Zanga

clear; clc; close all;

%% Settings
sourceToDetector = 119.40;    % mm
sourceToCenter = 105.80;      % mm
pixelSpacing = 0.020;        % mm per detector pixel
imageSpacing = 0.080;        % mm per reconstructed pixel
imageWidth = 350;            % reconstructed image columns
imageHeight = 300;           % reconstructed image rows
rotationDirection = -1;      % scanner rotation convention
crystalRow = 750;            % row to reconstruct
insectRow = 0;               % 0 means use the middle image row
crystal180Rotation = 1;      % 1 = rotate 90 degrees counterclockwise

%% Exercise 1: Crystal 360
disp('Select Crystal 360 Folder')
folder360 = uigetdir(pwd);

[crystal360Sino, crystal360Angles] = loadImages( ...
    folder360, '.bmp', crystalRow, 0);

showSino(crystal360Sino, crystal360Angles, 'Crystal 360 Sinogram');

crystal360Image = fbp(crystal360Sino, crystal360Angles, ...
    sourceToDetector, sourceToCenter, pixelSpacing, imageSpacing, ...
    imageWidth, imageHeight, rotationDirection);

showImage(crystal360Image, 'Crystal 360 Reconstruction');

%% Exercise 2: Insect Head 360
disp('Select Insect Head Folder')
folderInsect = uigetdir(pwd);

[insectSino, insectAngles] = loadImages( ...
    folderInsect, '.tif', insectRow, 0);
showSino(insectSino, insectAngles, 'Insect Head Sinogram');

insectImage = fbp(insectSino, insectAngles, ...
    sourceToDetector, sourceToCenter, pixelSpacing, imageSpacing, ...
    imageWidth, imageHeight, rotationDirection);

showImage(insectImage, 'Insect Head Reconstruction');

%% Exercise 3: Crystal 180
disp('Select Crystal 180 Folder')
folderCrystal180 = uigetdir(pwd);

[crystal180Sino, crystal180Angles] = loadImages( ...
    folderCrystal180, '.bmp', crystalRow, crystal180Rotation);

showSino(crystal180Sino, crystal180Angles, 'Crystal 180 Sinogram');

crystal180Image = fbp(crystal180Sino, crystal180Angles, ...
    sourceToDetector, sourceToCenter, pixelSpacing, imageSpacing, ...
    imageWidth, imageHeight, rotationDirection);

showImage(crystal180Image, 'Crystal 180 Reconstruction');

%% Function 1: Load projection images and create a sinogram
function [sinogram, angles] = loadImages(folder, extension, sliceRow, turns)

    files = dir(fullfile(folder, ['*' extension]));
    if isempty(files)
        error('No %s images found in folder: %s', extension, folder);
    end

    % Read the projection angle from each filename.
    % Example: Crystal_24.bmp corresponds to 24 degrees.
    angles = zeros(length(files), 1);
    for k = 1:length(files)
        filename = files(k).name;
        part = regexp(filename, '([0-9]+\.?[0-9]*)\.[^.]+$', 'tokens', 'once');
        if isempty(part)
            error('Cannot read angle from filename: %s', filename);
        end
        angles(k) = str2double(part{1});
    end

    % Sort projections into increasing angle order.
    [angles, order] = sort(angles);
    files = files(order);

    % Remove a duplicate image acquired at exactly 360 degrees.
    if angles(end) == 360
        angles(end) = [];
        files(end) = [];
    end

    first = imread(fullfile(folder, files(1).name));
    if size(first, 3) == 3
        first = rgb2gray(first);
    end
    first = rot90(first, turns);
    [imageRows, detectorCount] = size(first);

    if sliceRow == 0
        sliceRow = round(imageRows / 2);
    end
    if sliceRow < 1 || sliceRow > imageRows
        error('Slice row %d is outside image height %d.', sliceRow, imageRows);
    end

    sinogram = zeros(detectorCount, length(files));
    for k = 1:length(files)
        img = imread(fullfile(folder, files(k).name));
        if size(img, 3) == 3
            img = rgb2gray(img);
        end
        img = rot90(img, turns);
        if size(img, 1) ~= imageRows || size(img, 2) ~= detectorCount
            error('Image sizes do not match: %s', files(k).name);
        end
        sinogram(:, k) = double(img(sliceRow, :))';
    end

    fprintf('Loaded %d images from %s, using row %d.\n', ...
        length(files), folder, sliceRow);
end

%% Function 2: Manual fan-beam filtered back projection
function image = fbp(sinogram, angles, sourceToDetector, ...
    sourceToCenter, pixelSpacing, imageSpacing, ...
    imageWidth, imageHeight, rotationDirection)

    [detectorCount, viewCount] = size(sinogram);
    if length(angles) ~= viewCount || viewCount < 2
        error('Sinogram and angles do not match.');
    end

    % Physical detector coordinates (mm).
    detector = ((1:detectorCount)' - (detectorCount+1)/2) * pixelSpacing;

    % Fan-beam weighting before filtering.
    weight = sourceToDetector ./ sqrt(sourceToDetector^2 + detector.^2);
    weightedSino = sinogram .* weight;

    % Build the Ram-Lak filter in the frequency domain.
    fftLength = 2^nextpow2(2*detectorCount);
    n = (-fftLength/2:fftLength/2-1)';
    ramp = zeros(fftLength, 1);
    ramp(fftLength/2+1) = 1/4;
    for j = 1:fftLength
        if mod(abs(n(j)), 2) == 1
            ramp(j) = -1 / (pi*n(j))^2;
        end
    end
    filter = real(fft(ifftshift(ramp))) * 2;
    filter = max(filter, 0);

    % Filter each projection separately.
    startIndex = fftLength/2 - floor(detectorCount/2) + 1;
    filteredSino = zeros(detectorCount, viewCount);
    for k = 1:viewCount
        padded = zeros(fftLength, 1);
        padded(startIndex:startIndex+detectorCount-1) = weightedSino(:, k);
        result = real(ifft(fft(padded) .* filter));
        filteredSino(:, k) = result(startIndex:startIndex+detectorCount-1) ...
            / (2*pixelSpacing) * (sourceToDetector/sourceToCenter);
    end

    % Create the reconstruction grid (mm).
    x = ((1:imageWidth) - (imageWidth+1)/2) * imageSpacing;
    y = ((1:imageHeight) - (imageHeight+1)/2) * imageSpacing;
    [X, Y] = meshgrid(x, y);
    image = zeros(imageHeight, imageWidth);

    % Angular weight for equally spaced projections.
    angleStep = angles(2) - angles(1);
    angleWeight = (angleStep*pi/180)/2;

    % Back-project filtered projections.
    for k = 1:viewCount
        theta = angles(k) * rotationDirection * pi/180;
        rotatedX = X*cos(theta-pi/2) + Y*sin(theta-pi/2);
        rotatedY = -X*sin(theta-pi/2) + Y*cos(theta-pi/2);

        distance = rotatedY + sourceToCenter;
        detectorPoint = sourceToDetector * rotatedX ./ distance;
        backWeight = (sourceToCenter ./ distance).^2;

        values = interp1(detector, filteredSino(:, k), ...
            -detectorPoint, 'linear', 0);
        image = image + backWeight .* values * angleWeight;
    end
end

%% Function 3: Display sinogram
function showSino(sinogram, angles, name)
    figure;
    imagesc(angles, 1:size(sinogram,1), sinogram);
    colormap gray;
    colorbar;
    xlabel('Projection angle (degrees)');
    ylabel('Detector position (pixels)');
    title(name);
end

%% Function 4: Display reconstructed image
function showImage(image, name)
    figure;
    imagesc(image);
    axis image;
    colormap gray;
    colorbar;
    xlabel('X position (pixels)');
    ylabel('Y position (pixels)');
    title(name);
end
