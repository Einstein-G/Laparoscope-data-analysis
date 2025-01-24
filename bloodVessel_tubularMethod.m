close all;
clc;
clear;
warning('off');
names={'1749','2616','2801','4390','4279','1450','4928','4978','3425','4155','4551','1341','1311','3936','3889','3109','2972','2891','2645','2347','2283','0023','0037','2615','2769','3321','3388','2265','3631','0064','0503','0432','1413','1848'};
j=0;
fprintf('\t\tLesion\t\t\t\tHeatlthy\n');

set(0,'DefaultFigureWindowStyle','normal')

for i = 2:2%length(names)%patient2:1-35, patient4:36- ,patient6:


filename=names{i};
matfilename=sprintf('%s%s%s','t0',names{i},'.mat');
load(matfilename);
filename1 = sprintf('%s%s%s','diff_t0',filename,'_c00001.tif'); %diff red
filename2 = sprintf('%s%s%s','diff_t0',filename,'_c00002.tif'); %diff green
filename3 = sprintf('%s%s%s','diff_t0',filename,'_c00003.tif'); %diff blue
filename4 = sprintf('%s%s%s','sum_t0',filename,'_c00001.tif'); %sum red
filename5 = sprintf('%s%s%s','sum_t0',filename,'_c00002.tif'); %sum green
filename6 = sprintf('%s%s%s','sum_t0',filename,'_c00003.tif'); %sum blue

rd = imread(filename1); % red diff
gd = imread(filename2); 
bd = imread(filename3);
rs = imread(filename4); % red sum
gs = imread(filename5);
bs = imread(filename6);
rp= (rs-rd)./2; % red perpendicular
gp= (gs-gd)./2;
bp= (bs-bd)./2;

xi=51;xm=size(rd,1)-50;yi=51;ym=size(rd,2)-50;
%%%%%%%%glare mask

diff = cat(3,rd,gd,bd);
sum = cat(3,rs,gs,bs);
perp = cat(3,rp, gp, bp);

sigma = 4;
perpSum = imgaussfilt(normalise(double(rp + gp + bp)),sigma);
diffSum = imgaussfilt(normalise(double(rd + gd + bd)),sigma);

% CHOOSE IMAGE TYPE HERE

grayDiff = normalise(double(rgb2gray(diff)));
grayPerp = normalise(double(rgb2gray(perp)));

%% Manually selected thresholds
perpThresh = [0.05, 0.26];
diffThresh = [0.0500, 0.192];

[perpSumMask, perpSumEdges] = glareMask(perpSum, perpThresh, 5);
[diffSumMask, diffSumEdges] = glareMask(diffSum, diffThresh, 5);

overlapEdges = perpSumEdges + diffSumEdges - 1; % merging binary images
overlapMask = perpSumMask + diffSumMask - 1;
overlapMask(overlapMask < 0) = 0; % saturate negative values
overlapMask = normalise(overlapMask);

% Identify dark areas in overlapMask after burning onto grayPerp (lesions
% are dark.
getBlack = grayPerp;
darkThresh = prctile(grayPerp, 0.5,'all');
getBlack(~overlapMask) = 1;
getBlack(grayPerp >= darkThresh) = 1;
getBlack = imbinarize(getBlack);
getBlack = ~getBlack;
%% 
% Identify bright areas in overlapMask after burning onto grayDiff (lesions
% are bright.
getBright = grayDiff;
brighThresh = prctile(grayDiff, 99.5,'all');
getBright(~overlapMask) = 0;
getBright(grayDiff <= brighThresh) = 0;
getBright = imbinarize(getBright);

% if at least either of getbright or getdark are empty, return nothing
% as we expect to have both dark and white spots on glare.
% Assumption is made from observations, lesion are black in perp and bright
% in diff
g = zeros(size(grayPerp));
if (length(getBlack(getBlack==1)) > 50) & (length(getBright(getBright==1)) >= 50)
    g = imbinarize(getBlack + getBright);
    g = imdilate(g,strel('disk',2));
    g = imfill(g,'holes'); 
%     glare = bwconvhull(glare,'objects');
    g = bwareaopen(g,50);
end

%%%%%%%%
rd = normalise(im2double(rd)); % red diff
gd = normalise(im2double(gd));
bd = normalise(im2double(bd));
rs = normalise(im2double(rs)); % red sum
gs = normalise(im2double(gs));
bs = normalise(im2double(bs));
rp=  normalise(im2double(rp)); % red perpendicular
gp=  normalise(im2double(gp));
bp=  normalise(im2double(bp));

rd2 = glare(rd); % red diff
gd2 = glare(gd);
bd2 = glare(bd);
rs2 = glare(rs); % red sum
gs2 = glare(gs);
bs2 = glare(bs);
rp2 = glare(rp); % red perpendicular
gp2 = glare(gp);
bp2 = glare(bp);

bp3=bp2(51:size(bp2,1)-50,51:size(bp2,2)-50);
gp3=gp2(51:size(gp2,1)-50,51:size(gp2,2)-50);
rp3=rp2(51:size(rp2,1)-50,51:size(rp2,2)-50);
%% Blood Vessel Segmentation based on Green Channel
figure,montage({rp,gp,bp}), title("{rp,gp,bp}");
figure,montage({rescale(rp3),rescale(gp3),rescale(bp3)}), title("{rp3,gp3,bp3}");


%% Applying Fibermetrix filter with illumination corrected channels
red = fibermetricFilter(rp3, 5); % only extract 5 largest blobs
green = fibermetricFilter(gp3, 5);
blue = fibermetricFilter(bp3, 5);
figure, montage({red, green,blue}), title("segmented RGB mask") 

result = green; % Only using green and blue channels
result(green < 1) = 0;
result(blue < 1) = 0;
figure, imshow(result), title('result')

end

function output = fibermetricFilter(I, numBlob)
% this function takes input image and applies fibermitric filter to extract 
% lines and curves of width 40,50,60,70, binarizing and extracting numBlob 
% number of largest blobs
B = rescale(fibermetric(I,[40 50 60 70],'ObjectPolarity','dark'));
C= imbinarize(rescale(B), 0.1);
output = bwareafilt(C, numBlob); 

end
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

function output = normalise(image)
%normalize normalizes image matrix between 0 and 1
v = image - min(image(:));
v = v ./ max(v(:));
output = v;
end

function output = glare(image)
%illuminationCorrect correct illumination with fourier transformation
% correcting nonuniform ilumination, reference: https://blogs.mathworks.com/steve/2013/06/25/homomorphic-filtering-part-1/
I = log(1 + double(image));
M = 2*size(I,1) + 1;
N = 2*size(I,2) + 1;
sigma = 15;
[X, Y] = meshgrid(1:N,1:M);
centerX = ceil(N/2);
centerY = ceil(M/2);
gaussianNumerator = (X - centerX).^2 + (Y - centerY).^2;
H = exp(-gaussianNumerator./(2*sigma.^2));
H = 1 - H;
H = fftshift(H);
If = fft2(I, M, N);
Iout = real(ifft2(H.*If));
Iout = Iout(1:size(I,1),1:size(I,2));
output = exp(Iout) - 1;
end

% https://www.mathworks.com/matlabcentral/fileexchange/24990-retinal-blood-vessel-extraction
% https://www.researchgate.net/publication/334096859_Blood_Vessel_Extraction_Using_Combination_of_Kirsch's_Templates_and_Fuzzy_C-Means_FCM_on_Retinal_Images
function bloodVessels = VesselExtract(inImg, threshold)
                                                                           %Kirsch's Templates
h1=[5 -3 -3;
    5  0 -3;
    5 -3 -3]/15;
h2=[-3 -3 5;
    -3  0 5;
    -3 -3 5]/15;
h3=[-3 -3 -3;
     5  0 -3;
     5  5 -3]/15;
h4=[-3  5  5;
    -3  0  5;
    -3 -3 -3]/15;
h5=[-3 -3 -3;
    -3  0 -3;
     5  5  5]/15;
h6=[ 5  5  5;
    -3  0 -3;
    -3 -3 -3]/15;
h7=[-3 -3 -3;
    -3  0  5;
    -3  5  5]/15;
h8=[ 5  5 -3;
     5  0 -3;
    -3 -3 -3]/15;
%Spatial Filtering by Kirsch's Templates
t1=filter2(h1,inImg);
t2=filter2(h2,inImg);
t3=filter2(h3,inImg);
t4=filter2(h4,inImg);
t5=filter2(h5,inImg);
t6=filter2(h6,inImg);
t7=filter2(h7,inImg);
t8=filter2(h8,inImg);
s=size(inImg);
bloodVessels=zeros(s(1),s(2));
temp=zeros(1,8);
for i=1:s(1)
    for j=1:s(2)
        temp(1)=t1(i,j);temp(2)=t2(i,j);temp(3)=t3(i,j);temp(4)=t4(i,j);
        temp(5)=t5(i,j);temp(6)=t6(i,j);temp(7)=t7(i,j);temp(8)=t8(i,j);
        if(max(temp)>threshold)
            bloodVessels(i,j)=max(temp);
        end
    end
end
end