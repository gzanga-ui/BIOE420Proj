%% BIOE420 Project 1: CT Image Reconstruction 
% Group Members: Julie Karasik, Gabrielle Zanga

clear; clc; close all;


scanner.DSD = 119.40;     % mm, source to detector
scanner.DSO = 105.80;     % mm, source to rotation axis
scanner.du  = 0.020;      % mm per detector pixel
scanner.dir = -1;         % instructor's rotation direction
scanner.nx  = 350;        % reconstructed image width
scanner.ny  = 300;        % reconstructed image height
scanner.dx  = 0.080;      % mm per reconstruction pixel
scanner.dy  = 0.080;      % mm per reconstruction pixel
scanner.filter = 'ram-lak';

% If detector center is offset from the ideal midpoint, calibrate it here.
% With 1000 detector pixels, geometric center is 500.5 (between pixels).
% The assignment's (500,750) identifies the approximate slice location.

%% =========================================================
% EXERCISE 1: CRYSTAL DATASET, 360 DEGREES
% ==========================================================
[crystal360Sinogram, crystal360Angles] = loadProjectionLines( ...
    'Data_Crystal_360.zip', 'Crystal 360 Images', '.bmp', 750, 0);
figure;
showSinogram(crystal360Sinogram, crystal360Angles, ...
    'Exercise 1: Crystal 360-Degree Sinogram');
reconstructedCrystal360 = reconstructFanBeam2D( ...
    crystal360Sinogram, crystal360Angles, scanner);
figure;
showReconstruction(reconstructedCrystal360, ...
    'Exercise 1: Crystal 360-Degree Fan-Beam FBP');

%% =========================================================
% EXERCISE 2: INSECT HEAD DATASET, 360 DEGREES
% ==========================================================
% For insect images, [] chooses the middle row automatically.
[insectSinogram, insectAngles] = loadProjectionLines( ...
    'Dataset_Insect_Head_360.zip', 'Insect Head Images', '.tif', [], 0);
figure;
showSinogram(insectSinogram, insectAngles, ...
    'Exercise 2: Insect Head 360-Degree Sinogram');
reconstructedInsect = reconstructFanBeam2D( ...
    insectSinogram, insectAngles, scanner);
figure;
showReconstruction(reconstructedInsect, ...
    'Exercise 2: Insect Head Fan-Beam FBP');

%% =========================================================
% EXERCISE 3: CRYSTAL DATASET, 180 DEGREES
% ==========================================================
% Instructor requires rotating the raw images 90 degrees.
% +1 = 90 degrees counterclockwise; set -1 if raw image requires clockwise.
crystal180Rotation = 1;
[crystal180Sinogram, crystal180Angles] = loadProjectionLines( ...
    'Dataset_Crystal_180.zip', 'Crystal_180_Images', '.bmp', ...
    750, crystal180Rotation);
figure;
showSinogram(crystal180Sinogram, crystal180Angles, ...
    'Exercise 3: Crystal 180-Degree Sinogram');
reconstructedCrystal180 = reconstructFanBeam2D( ...
    crystal180Sinogram, crystal180Angles, scanner);
figure;
showReconstruction(reconstructedCrystal180, ...
    'Exercise 3: Crystal 180-Degree Fan-Beam FBP (Limited Angle)');


%% =========================================================
% LOCAL FUNCTIONS (MUST BE AT END OF THIS SCRIPT)
% ==========================================================
function [sino, angles] = loadProjectionLines(zipName, folder, extension, sliceRow, rotateCount)
    if ~isfolder(folder)
        if ~isfile(zipName)
            error('Cannot find %s. Put the ZIP in MATLAB Current Folder.',zipName);
        end
        unzip(zipName,folder);
    end
    if strcmpi(extension,'.tif')
        files = [dir(fullfile(folder,'**','*.tif')); ...
                 dir(fullfile(folder,'**','*.tiff'))];
    else
        files = dir(fullfile(folder,'**',['*' extension]));
    end
    if isempty(files)
        error('No projection images found under %s.',folder);
    end

    % Read angle immediately before extension, ignoring optional - Copy.
    anglesAll = nan(numel(files),1);
    for k=1:numel(files)
        tok = regexp(files(k).name, ...
            '([+-]?\d+(?:\.\d+)?)\s*(?:-\s*Copy)?\.(?:bmp|tiff?)$', ...
            'tokens','once','ignorecase');
        if ~isempty(tok)
            anglesAll(k) = mod(str2double(tok{1}),360);
        end
    end
    keep = isfinite(anglesAll);
    files = files(keep); anglesAll = anglesAll(keep);
    if isempty(files), error('Could not parse projection angles in %s',folder); end

    % Prefer original files to Windows-generated "- Copy" files.
    isCopy = contains({files.name},'- Copy','IgnoreCase',true);
    [~,ord] = sortrows([anglesAll, double(isCopy(:))],[1 2]);
    files = files(ord); anglesAll = anglesAll(ord);
    [angles,idx] = unique(anglesAll,'stable');
    files = files(idx);
    angles = angles(:);

    first = imread(fullfile(files(1).folder,files(1).name));
    if size(first,3)==3, first=rgb2gray(first); end
    first = rot90(first,rotateCount);
    [nRows,nDetectors] = size(first);
    if isempty(sliceRow), sliceRow=round(nRows/2); end
    if sliceRow<1 || sliceRow>nRows
        error('Slice row %d is outside the rotated image (height %d).',sliceRow,nRows);
    end
    sino = zeros(nDetectors,numel(angles));
    for k=1:numel(angles)
        im = imread(fullfile(files(k).folder,files(k).name));
        if size(im,3)==3, im=rgb2gray(im); end
        im = rot90(im,rotateCount);
        if size(im,1)~=nRows || size(im,2)~=nDetectors
            error('Inconsistent image size: %s',files(k).name);
        end
        sino(:,k)=double(im(sliceRow,:)).';
    end
    fprintf('%s: %d views, %g to %g deg, %d detectors, row %d, rot90=%d\n', ...
        folder,numel(angles),min(angles),max(angles),nDetectors,sliceRow,rotateCount);
    if any(diff(angles)>2.1)
        warning('%s has angle gaps larger than 2 degrees.',folder);
    end
end

function recon = reconstructFanBeam2D(sino, angles, cfg)
    % Fan-beam FBP central-slice adaptation of the instructor's FDK code.
    % Assumes line-integral projections and a FLAT detector.
    % If image values are raw intensities, -log(I/I0) preprocessing may be
    % required before reconstruction; I0 calibration is not provided here.
    [nu,nViews] = size(sino);
    angles = angles(:);
    if numel(angles)~=nViews, error('Angles / sinogram mismatch.'); end
    if nViews<2, error('At least two projections are needed.'); end

    du=cfg.du; DSD=cfg.DSD; DSO=cfg.DSO;
    % Physical detector positions (mm) about panel midpoint.
    detectorU = ((1:nu)'-(nu+1)/2)*du;
    % Instructor's detector pre-weighting, evaluated at central slice v=0.
    preweight = DSD./sqrt(DSD^2+detectorU.^2);
    weighted = sino.*repmat(preweight,1,nViews);

    % Professor's spatial-domain ramp kernel -> FFT frequency response.
    % Padding and zero-centering avoid wraparound at the detector edges.
    L = 2^nextpow2(2*nu);
    nn = (-L/2:L/2-1)';
    kernel = zeros(L,1);
    kernel(L/2+1) = 1/4;
    odd = mod(abs(nn),2)==1;
    kernel(odd) = -1./(pi*nn(odd)).^2;
    H = real(fft(ifftshift(kernel)))*2;
    H = max(H,0);  % Ram-Lak, no extra Hann/Hamming apodization
    padStart = L/2-floor(nu/2)+1;
    filtered = zeros(nu,nViews);
    for k=1:nViews
        padded = zeros(L,1);
        padded(padStart:padStart+nu-1)=weighted(:,k);
        fp = real(ifft(fft(padded).*H));
        % Same detector-spacing and geometric scaling form as instructor.
        filtered(:,k)=fp(padStart:padStart+nu-1) ...
            /(2*du)*(DSD/DSO);
    end

    % Voxel grid (mm), mirroring ParamSettingr16F.m.
    xs=((1:cfg.nx)-(cfg.nx+1)/2)*cfg.dx;
    ys=((1:cfg.ny)-(cfg.ny+1)/2)*cfg.dy;
    [xx,yy]=meshgrid(xs,ys);
    recon=zeros(cfg.ny,cfg.nx);

    % Angular integration. Actual 180-degree fan-beam scan is incomplete:
    % no missing angles are synthesized, and Parker weighting is NOT used.
    if max(angles)-min(angles)>270
        % Circular quadrature including wraparound for complete 360 scan.
        prev=[angles(end)-360;angles(1:end-1)];
        next=[angles(2:end);angles(1)+360];
        weights=deg2rad((next-prev)/2)/2;
    else
        % Limited-angle integral over the acquired range (2-degree nominal).
        step=median(diff(angles));
        left=[angles(1)-step/2; (angles(1:end-1)+angles(2:end))/2];
        right=[(angles(1:end-1)+angles(2:end))/2; angles(end)+step/2];
        weights=deg2rad(right-left)/2;
        warning(['This 180-degree fan-beam acquisition is an incomplete ', ...
                 'short scan. Strong streak artifacts are expected.']);
    end

    for k=1:nViews
        % Match instructor's backprojection_line.m coordinate convention.
        a=deg2rad(angles(k)*cfg.dir);
        rx=xx*cos(a-pi/2)+yy*sin(a-pi/2);
        ry=-xx*sin(a-pi/2)+yy*cos(a-pi/2);
        denominator=ry+DSO;
        u=DSD*rx./denominator;
        weight=(DSO./denominator).^2;
        % Professor maps projected coordinate to detector index with a
        % negative detector axis. Use interpolation with zero outside FOV.
        sample=interp1(detectorU,filtered(:,k),-u,'linear',0);
        sample(~isfinite(sample))=0;
        recon=recon+weight.*sample.*weights(k);
    end
end

function showSinogram(sino,angles,ttl)
    imagesc(angles,1:size(sino,1),sino);
    colormap(gray); colorbar; axis normal;
    xlabel('Projection Angle (degrees)');
    ylabel('Detector Position (pixels)'); title(ttl);
end

function showReconstruction(img,ttl)
    imagesc(img); colormap(gray); axis image; colorbar;
    xlabel('X Position (pixels)'); ylabel('Y Position (pixels)');
    title(ttl);
end
