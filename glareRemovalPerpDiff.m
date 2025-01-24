clc;
clear;
clf; close;
clc;  % Clear command window.
clear;  % Delete all variables.
close all;  % Close all figure windows except those created by imtool.
imtool close all;  % Close all figure windows created by imtool.
workspace;  % Make sure the workspace panel is showing.
names={'1749','2616','2801','4390','4279','4279','1450','2801','4928','4978','4978','4978','3425','3425','4155','4551','1341','1341','1311','3936','3889','3109','2972','2972','2891','2891','2891','2891','2891','2645','2645','2347','2347','2283','2283'};
names2={'1749','2616','2801','4390','4279','4279_2','1450','2801_2','4928','4978_2','4978_3','4978','3425','3425_2','4155','4551','1341','1341_2','1311','3936','3889','3109','2972','2972_2','2891','2891_2','2891_3','2891_4','2891_5','2645','2645_2','2347','2347_2','2283','2283_2'};
fontSize = 12;
% % Ask user for image path input
prompt = {'Select image number:'};
dlgtitle = 'Specify Input';
dims = [1 35];
definput = {'1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
i = str2double(fullfile(answer{1,1}));


filename=names{i};
% matfilename=sprintf('%s%s%s','t0',names2{i},'.mat');
% load(matfilename);
filename1 = sprintf('%s%s%s','diff_t0',filename,'_c00001.tif'); %diff red
filename2 = sprintf('%s%s%s','diff_t0',filename,'_c00002.tif'); %diff green
filename3 = sprintf('%s%s%s','diff_t0',filename,'_c00003.tif'); %diff blue
filename4 = sprintf('%s%s%s','sum_t0',filename,'_c00001.tif'); %sum red
filename5 = sprintf('%s%s%s','sum_t0',filename,'_c00002.tif'); %sum green
filename6 = sprintf('%s%s%s','sum_t0',filename,'_c00003.tif'); %sum blue

% Reading each image
rd = imread(filename1); % red diff
gd = imread(filename2); 
bd = imread(filename3);
rs = imread(filename4); % red sum
gs = imread(filename5);
bs = imread(filename6);
rp= (rs-rd)./2; % red perpendicular
gp= (gs-gd)./2;
bp= (bs-bd)./2;


diff = cat(3,rd,gd,bd);
sum = cat(3,rs,gs,bs);
perp = cat(3,rp, gp, bp);

sigma = 4;
perpSum = imgaussfilt(normalize(double(rp + gp + bp)),sigma);
diffSum = imgaussfilt(normalize(double(rd + gd + bd)),sigma);

% CHOOSE IMAGE TYPE HERE


grayDiff = normalize(double(rgb2gray(diff)));
grayPerp = normalize(double(rgb2gray(perp)));
% grayPerpMinusDiff = normalize(double(illuminationCorrect(double(rgb2gray(perp - diff)))));
% grayDiffMinusPerp = normalize(double(illuminationCorrect(double(rgb2gray(diff - perp)))));

figure
imshow(grayDiff),title('grayDiff'), axis on
figure
imshow(grayPerp),title('grayPerp'), axis on


%% Manually selected thresholds
perpThresh = [0.05, 0.26];
diffThresh = [0.0500, 0.192];

[perpSumMask, perpSumEdges] = glareMask(perpSum, perpThresh, 5);
[diffSumMask, diffSumEdges] = glareMask(diffSum, diffThresh, 5);

overlapEdges = perpSumEdges + diffSumEdges - 1; % merging binary images
overlapMask = perpSumMask + diffSumMask - 1;
overlapMask(overlapMask < 0) = 0; % saturate negative values
overlapMask = normalize(overlapMask);


figure,montage({perpSumMask, diffSumMask,overlapMask, diff});
title('mask');

% Identify dark areas in overlapMask after burning onto grayPerp (lesions
% are dark.
getBlack = grayPerp;
darkThresh = prctile(grayPerp, 0.5,'all');
getBlack(~overlapMask) = 1;
getBlack(grayPerp >= darkThresh) = 1;
getBlack = imbinarize(getBlack);
getBlack = ~getBlack;
figure, imshow(getBlack); title('get black lesion in Perp');
%% 
% Identify bright areas in overlapMask after burning onto grayDiff (lesions
% are bright.
getBright = grayDiff;
brighThresh = prctile(grayDiff, 99.5,'all');
getBright(~overlapMask) = 0;
getBright(grayDiff <= brighThresh) = 0;
getBright = imbinarize(getBright);
figure, imshow(getBright);title('get bright');

% if at least either of getbright or getdark are empty, return nothing
% as we expect to have both dark and white spots on glare.
% Assumption is made from observations, lesion are black in perp and bright
% in diff
glare = zeros(size(grayPerp));
if (length(getBlack(getBlack==1)) > 50) & (length(getBright(getBright==1)) >= 50)
    glare = imbinarize(getBlack + getBright);
    glare = imdilate(glare,strel('disk',2));
    glare = imfill(glare,'holes'); 
%     glare = bwconvhull(glare,'objects');
    glare = bwareaopen(glare,50);
  
end

figure('Name','Montage of Gray Diff, Perp, Glare Mask'),
multi = cat(3,grayDiff,grayPerp,glare);
montage(multi);

%GLARE MASK FUNCTION
%This function serves to apply Canny edge detection to identify glare in
%the supplied grayscale image and threshold value, then create a mask which
%describes the detected glare regions.
function [mask, fatEdges] = glareMask(image, threshold,dilation)
[BW, threshOut] = edge(image,'canny',threshold);
fatEdges = imdilate(BW,strel('disk',dilation));
% adding white edges for imfill
fatEdges(:,1) = 1;
fatEdges(:,size(fatEdges,2)) = 1;
mask = imfill(fatEdges,'holes'); %Mask returned by the function which includes
mask(:,1) = 0;
mask(:,size(mask,2)) = 0;
end


function output = normalize(image)
%normalize normalizes image matrix between 0 and 1
v = image - min(image(:));
v = v ./ max(v(:));
output = v;
end