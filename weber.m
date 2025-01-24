function [L3] = weber(w)
warning('off');
names={'1749','2616','2801','4390','4279','1450','4928','4978','3425','4155','4551','1341','1311','3936','3889','3109','2972','2891','2645','2347','2283','0023','0037','2615','2769','3321','3388','2265','3631','0064','0503','0432','1413','1848'};
j=0;
for i = 1:1%length(names)%patient2:1-35, patient4:36- ,patient6:


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

% rd2 = gaussianbpf_dp(rd2,25,256);
% rd2 = gaussianbpf_dp(rd2,20,120);
% rd2 = butterworthbpf_caroline(rd2,10,120,3);
% 
% gd2 = gaussianbpf_dp(gd2,25,256);
% gd2 = gaussianbpf_dp(gd2,20,120);
% gd2 = butterworthbpf_caroline(gd2,10,120,3);
% 
% bd2 = gaussianbpf_dp(bd2,25,256);
% bd2 = gaussianbpf_dp(bd2,20,120);
% bd2 = butterworthbpf_caroline(bd2,10,120,3);
% 
% rp2 = gaussianbpf_dp(rp2,25,256);
% rp2 = gaussianbpf_dp(rp2,20,120);
% rp2 = butterworthbpf_caroline(rp2,10,120,3);
% 
% gp2 = gaussianbpf_dp(gp2,25,256);
% gp2 = gaussianbpf_dp(gp2,20,120);
% gp2 = butterworthbpf_caroline(gp2,10,120,3);
% 
% bp2 = gaussianbpf_dp(bp2,25,256);
% bp2 = gaussianbpf_dp(bp2,20,120);
% bp2 = butterworthbpf_caroline(bp2,10,120,3);
% 
% rs2 = gaussianbpf_dp(rs2,25,256);
% rs2 = gaussianbpf_dp(rs2,20,120);
% rs2 = butterworthbpf_caroline(rs2,10,120,3);
% 
% gs2 = gaussianbpf_dp(gs2,25,256);
% gs2 = gaussianbpf_dp(gs2,20,120);
% gs2 = butterworthbpf_caroline(gs2,10,120,3);
% 
% bs2 = gaussianbpf_dp(bs2,25,256);
% bs2 = gaussianbpf_dp(bs2,20,120);
% bs2 = butterworthbpf_caroline(bs2,10,120,3);

var11=normalise(double(rp2-w*bp2));
var1=normalise((imcomplement(var11)));
var21=normalise(double((rd2-w*bd2)));
var2=normalise(imcomplement(var21));

var3=normalise(double(imadjust((rp2./bp2))));
var41=normalise(double(rd2./bd2));
var4=normalise(double(imcomplement(var41)));

%vart(~bwL4)=0;

backgroundImage = normalise(imadjust(rd+gd+bd));
backgroundImagec = backgroundImage(xi:xm,yi:ym);

% w     = 50;       % bilateral filter half-width5
% sigma = [50 0.5]; % bilateral filter standard deviations3 0.1

% Apply bilateral filter to each image.
% rdt = bfilter2(rd,w,sigma); % red diff

heatmapImage1 = imgaussfilt(var11,1);
heatmapImage2 = imgaussfilt(var21,5);
heatmapImage3 = imgaussfilt(var3,5);
heatmapImage4 = imgaussfilt(var41,5);
heatmapImage5=imadjust(im2double(0.2*heatmapImage1+0.4*heatmapImage2+0.8*heatmapImage3));

heatmapImage=normalise(imadjust(heatmapImage2));
g=imgaussfilt(heatmapImage,30).*g;
heatmapImage(g~=0)=g(g~=0);
heatmapImagec= heatmapImage(xi:xm,yi:ym);

sizeHeatmapImage = size(heatmapImagec);
sizeBackgroundImage = size(backgroundImagec);

heatmapTransparency = ones(size(backgroundImagec,1),size(backgroundImagec,2));
heatmapTransparency = heatmapTransparency*1;
%heatmapTransparency(~bwL4)=0;

bw1=bwL1(xi:xm,yi:ym);
bound = bwboundaries(bw1,'noholes');
n=size(bound,1);
hold on;
for k=1:n
    thisbound=bound{k};
    x=thisbound(:,2);
    y=thisbound(:,1);
    polyshapeL1 = polyshape(x,y);
areaL1 = area(polyshapeL1);
[x1,y1] = boundary(polyshapeL1);
[cx1, cy1] = centroid(polyshapeL1); %find the centroid

n = 1;
polyshapeL2 = 0;
areaL2 = 0;
while (areaL2 - areaL1) < (1.1*areaL1)
    polyshapeL2 = polybuffer(polyshapeL1, n); %expand by n units
    areaL2 = area(polyshapeL2);
    n = n+0.1;
end
[x2,y2] = boundary(polyshapeL2);

% Calculating L3
n = 1;
polyshapeL3 = 0;
areaL3 = 0;
while (areaL3 - areaL2) < (1*areaL1)
    polyshapeL3 = polybuffer(polyshapeL2, n); %expand by n units
    areaL3 = area(polyshapeL3);
    n = n+0.1;
end
[x3,y3] = boundary(polyshapeL3);

bwL1 = poly2mask(x1, y1, size(rs, 1), size(rs, 2)); %returns a binary mask size as input img
bwL2 = poly2mask(x2, y2, size(rs, 1), size(rs, 2)); %returns a binary mask size as input img
bwL3 = poly2mask(x3, y3, size(rs, 1), size(rs, 2)); %returns a binary mask size as input img
bwL4 = bwL3;
% bwL3(bwL2)=0;
    temp1=heatmapImage;temp3=heatmapImage;
    temp1(~bwL1)=0;
    temp3(~bwL3)=0;
    L1=mean2(nonzeros(temp1));
    L3=mean2(nonzeros(temp3));
end
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
end

