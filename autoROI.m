% Referrence code: https://www.mathworks.com/matlabcentral/answers/58000-how-can-i-draw-on-an-image-and-use-what-i-draw-as-segmentation-seeds
% Demo to have the user freehand draw an irregular shape over
% a gray scale image to indicate L1 (Region of interest ROI).
% Determine buffer boundary (L2) 10 pixels away from L1, and background
% boundary (L3) so that region enclosed by L3 and L2 has close area value
% to ROI area.
% NOTE: Comment line 42 if you are working with black white images to speed up.
% Luminance is obtained by averaging light intensity over region of
% interest.
% L3 is found by incrementing expansion unit by 0.1 from L2 while chbwL3eckingbwL3
% for background area (between L3 and L2) to be smaller than ROI. Will
% output
% the first area that is larger than ROI's area.

clc;  % Clear command window.
clear;  % Delete all variables.
close all;  % Close all figure windows except those created by imtool.
imtool close all;  % Close all figure windows created by imtool.
workspace;  % Make sure the workspace panel is showing.
fontSize = 12;

% Ask user for image path input
prompt = {'Enter folder path:','Enter image file name with extension (.tif, .png,..):'};
dlgtitle = 'Specify Input';
dims = [1 35];
definput = {'C:\Users\einst\Desktop\laparascope\Paper Instrumentation\Trout Masters Imaging Data\Human Peritoneum\Image Data\Laparo Setup\MetsHigh2\Sum','blueSum.png'};
answer = inputdlg(prompt,dlgtitle,dims,definput)
folder = fullfile(answer{1,1})
baseFileName = answer{2,1}


% Get the full filename, with path prepended.
fullFileName = fullfile(folder, baseFileName);
% Check if file exists.
if ~exist(fullFileName, 'file')
  % File doesn't exist -- didn't find it there.  Check the search path for it.
  fullFileName = baseFileName; % No path this time.
  if ~exist(fullFileName, 'file')
    % Still didn't find it.  Alert user.
    errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
    uiwait(warndlg(errorMessage));
    return;
  end
end

% Read image
grayImage = imread(fullFileName);
%grayImage = im2gray(grayImage); % TO HANDLE RGB IMAGES, uncomment if needed

% Drawing L1
imshow(grayImage, []);
axis off;
title('Original Grayscale Image', 'FontSize', fontSize);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
message = sprintf('Left click and hold to begin drawing region of L1 - ROI.\nSimply lift the mouse button to finish');
uiwait(msgbox(message));
hFH = imfreehand();
% Create a binary image L1 ("mask") from the ROI object.
binaryImageL1 = hFH.createMask();
xy = hFH.getPosition;
% Get coordinates of the boundary of the freehand drawn region L1.
structBoundaries = bwboundaries(binaryImageL1);
xy=structBoundaries{1}; % Get n by 2 array of x,y coordinates.
x = xy(:, 2); % Columns.
y = xy(:, 1); % Rows.
% END OF DRAWING L1


% L1 boundary properties
polyshapeL1 = polyshape(x,y); %create a polygon from x y coordinates
areaL1 = area(polyshapeL1);;
[x1,y1] = boundary(polyshapeL1);
[cx1, cy1] = centroid(polyshapeL1); %find the centroid

% Calculating L2
polyshapeL2 = polybuffer(polyshapeL1,+10); %expand by 10 units
areaL2 = area(polyshapeL2);
[x2,y2] = boundary(polyshapeL2);

% Calculating L3
% The background region between L3 and L2 needs to have area similar to
% ROI's area (areaL1).
% Using while loop to increment expansion unit until difference between
% areaL3 and areaL2 = areaL1
n = 1;
polyshapeL3 = 0;
areaL3 = 0;
while areaL3 - areaL2 < areaL1
    polyshapeL3 = polybuffer(polyshapeL2, n); %expand by n units
    areaL3 = area(polyshapeL3)
    n = n+0.1;
end
[x3,y3] = boundary(polyshapeL3);


% Generating black & white masks
bwL1 = poly2mask(x1, y1, size(grayImage, 1), size(grayImage, 2)); %returns a binary mask size as input img
bwL2 = poly2mask(x2, y2, size(grayImage, 1), size(grayImage, 2)); %returns a binary mask size as input img
bwL3 = poly2mask(x3, y3, size(grayImage, 1), size(grayImage, 2)); %returns a binary mask size as input img

% Generating image of masked ROI and background areas, calculate mean
% luminance

% ROI (inside L1)
blackMaskedImageL1 = grayImage;
blackMaskedImageL1(~bwL1) = 0;
luminanceL1 = mean2(nonzeros(blackMaskedImageL1));
blackMaskedImageL2 = grayImage;
blackMaskedImageL2(~bwL2) = 0;
luminanceL2 = mean2(nonzeros(blackMaskedImageL2));
% Luminance of region of interest
% Inside L3
blackMaskedImageL3 = grayImage;
blackMaskedImageL3(~bwL3) = 0; %region enclosed by L3 to be image
blackMaskedImageL3(bwL2) = 0; %region enclosed by L2 to be 0
luminanceL3 = mean2(nonzeros(blackMaskedImageL3)); % Note this is background region, not entire L3

% Weber contrast (lesion - background) / Background
contrast = (luminanceL1 - luminanceL3)/luminanceL3;

% Signal to background ratio Lesion / Background
ratio = luminanceL1 / luminanceL3;

% Areas
backgroundArea = areaL3 - areaL2;
ROIarea = areaL1;
areaDifference = backgroundArea - areaL1;

imshow(grayImage);hold on; plot(x3,y3, 'r');
drawnow; % Force it to draw immediately.
plot(x2,y2, 'b');
drawnow; % Force it to draw immediately.
plot(x,y, 'g');
drawnow; % Force it to draw immediately.

% % PLOTTING
% % Plotting boundaries
% subplot(1,1,1);
% imshow(grayImage)%, 'InitialMagnification', 800);
% axis on;
% title('L1 (green) & L2 (blue) & L3 (red) outlines', 'FontSize', fontSize);
% 
% subplot(1,1,1); % Plot L2 over grayImage with L1.
% hold on; % Don't blow away the image.
% % L1
% plot(x3,y3, 'r');
% drawnow; % Force it to draw immediately.
% plot(x2,y2, 'b');
% drawnow; % Force it to draw immediately.
% plot(x,y, 'g');
% drawnow; % Force it to draw immediately.
% 
% 
% figure
% % Display background and region of interest masked image
% subplot(2,1,1)
% imshow(blackMaskedImageL1, 'InitialMagnification', 800)
% title(sprintf('Region of interest (L1), mean luminance: %f', luminanceL1))
% 
% subplot(2,1,2)
% imshow(blackMaskedImageL3, 'InitialMagnification', 800)
% title(sprintf('Background (L3-L2), mean luminance: %f', luminanceL3))
% 



message = sprintf('Weber Contrast= %.s\nLesion/Background ratio= %f\nArea ROI (L1) = %.f\nArea background (L3-L2) = %.f\nAreaDifference = %.f\nMean luminance 10px buffer region = %.f', ...
contrast, ratio, ROIarea, backgroundArea, areaDifference, luminanceL2)
msgbox(message);