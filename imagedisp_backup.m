close all;
clc;
clear;
names={'1749','2616','2801','4390','4279','4279','1450','2801','4928','4978','4978','4978','3425','3425','4155','4551','1341','1341','1311','3936','3889','3109','2972','2972','2891','2891','2891','2891','2891','2645','2645','2347','2347','2283','2283'};
names2={'1749','2616','2801','4390','4279','4279_2','1450','2801_2','4928','4978_2','4978_3','4978','3425','3425_2','4155','4551','1341','1341_2','1311','3936','3889','3109','2972','2972_2','2891','2891_2','2891_3','2891_4','2891_5','2645','2645_2','2347','2347_2','2283','2283_2'};

for i = 1:1%length(names)%patient2:1-35, patient4:36- ,patient6:
    
filename=names{i};
matfilename=sprintf('%s%s%s','t0',names2{i},'.mat');
% matfilename2=sprintf('%s%s%s','t0',names2{i+1},'.mat');
% matfilename3=sprintf('%s%s%s','t0',names2{i+2},'.mat');
load(matfilename);
% load(matfilename2);
% load(matfilename3);
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

rd = im2double(rd); % red diff
gd = im2double(gd);
bd = im2double(bd);
rs = im2double(rs); % red sum
gs = im2double(gs);
bs = im2double(bs);
rp=  im2double(rp); % red perpendicular
gp=  im2double(gp);
bp=  im2double(bp);

% Set bilateral filter parameters.
w     = 5;       % bilateral filter half-width
sigma = [3 0.1]; % bilateral filter standard deviations

% Apply bilateral filter to each image.
% rdt = bfilter2(rd,w,sigma); % red diff
% gdt = bfilter2(gd,w,sigma);
% bdt = bfilter2(bd,w,sigma);


% rs2 = bfilter2(rs1,w,sigma); % red sum
% gs2 = bfilter2(gs1,w,sigma);
% bs2 = bfilter2(bs1,w,sigma);
% rp2 = bfilter2(rp1,w,sigma); % red perpendicular
% gp2 = bfilter2(gp1,w,sigma);
% bp2 = bfilter2(bp1,w,sigma); 

% rdt1 = im2double(glare(rdt)); % red diff
% gdt1 = im2double(glare(gdt));
% bdt1 = im2double(glare(bdt));



rd2 = glare(rd); % red diff
gd2 = glare(gd);
bd2 = glare(bd);
rs2 = glare(rs); % red sum
gs2 = glare(gs);
bs2 = glare(bs);
rp2 = glare(rp); % red perpendicular
gp2 = glare(gp);
bp2 = glare(bp);

% rd2=(rd2-min(min(rd2)))./(max(max(rd2))-min(min(rd2)));
% gd2=(rd2-min(min(rd2)))./(max(max(rd2))-min(min(rd2)));
% bd2=(rd2-min(min(rd2)))./(max(max(rd2))-min(min(rd2)));
% rs2=(rd2-min(min(rd2)))./(max(max(rd2))-min(min(rd2)));
% gs2=(rd2-min(min(rd2)))./(max(max(rd2))-min(min(rd2)));
% bs2=(rd2-min(min(rd2)))./(max(max(rd2))-min(min(rd2)));
% rp2=(rd2-min(min(rd2)))./(max(max(rd2))-min(min(rd2)));
% gp2=(rd2-min(min(rd2)))./(max(max(rd2))-min(min(rd2)));
% bp2=(rd2-min(min(rd2)))./(max(max(rd2))-min(min(rd2)));

structBoundaries = bwboundaries(bwL1);
xy=structBoundaries{1}; % Get n by 2 array of x,y coordinates.
x = xy(:, 2); % Columns.
y = xy(:, 1);
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
bwL3(bwL2)=0;

var11=im2double((gp2-bp2));
var12=(var11-min(min(var11)))./(max(max(var11))-min(min(var11)));
var12=(imcomplement(var12));
var1=(var12-min(min(var12)))./(max(max(var12))-min(min(var12)));
var21=im2double((gd2-bd2));
var22=(var21-min(min(var21)))./(max(max(var21))-min(min(var21)));
var22=(imcomplement(var22));
var2=(var22-min(min(var22)))./(max(max(var22))-min(min(var22)));

var3=im2double(imadjust((gp2./bp2)));
var41=im2double(gd2./bd2);
var4=im2double(imcomplement(var41));

%vart(~bwL4)=0;

backgroundImage = imadjust(rd+gd+bd);
backgroundImagec = backgroundImage(51:1000,51:1100);
% percentage1 = 0;
% percentage2 =1;
% low_bound = percentage1 / 100;
% high_bound = 1 - percentage2/100;
% K = imadjust(Im1,stretchlim(Im1,[low_bound high_bound]),[]);


% heatmapImage1 = medfilt2(var1,[5 5]);
% heatmapImage2 = medfilt2(var2,[5 5]);
% heatmapImage3 = medfilt2(var3,[5 5]);
% heatmapImage4 = medfilt2(var41,[5 5]);
heatmapImage1 = imgaussfilt(var1,20);
heatmapImage2 = imgaussfilt(var21,20);
heatmapImage3 = imgaussfilt(var3,20);
heatmapImage4 = imgaussfilt(var41,4);

heatmapImage5=imadjust(im2double(0.2*heatmapImage1+0.4*heatmapImage2+0.8*heatmapImage3));
%heatmapImage5=(im2double(rs+gs+bs));
% mask=heatmapImage5(:,:)<=1;
% heatmapImage=heatmapImage5.*mask;
%  heatmapImage=(heatmapImage5);
%heatmapImage=imadjust(im2double(rd+gd+bd));
heatmapImage=(heatmapImage5-min(min(heatmapImage5)))./(max(max(heatmapImage5))-min(min(heatmapImage5)));
heatmapImagec= heatmapImage5(51:1000,51:1100);

v1=heatmapImage;
v2=heatmapImage;
v1(~bwL1)=0;
v2(~bwL3)=0;
L1=mean2(nonzeros(v1));
L3=mean2(nonzeros(v2));
weber(i)=(L1-L3)./L3;
fun=1/weber;

%heatmapImage = bfilter2(var,20,[8 0.3]);
sizeHeatmapImage = size(heatmapImagec);
sizeBackgroundImage = size(backgroundImagec);

heatmapTransparency = ones(size(backgroundImagec,1),size(backgroundImagec,2));
heatmapTransparency = heatmapTransparency*1;
%heatmapTransparency(~bwL4)=0;

% Showing background image and overlaying heatmap
figure;
ax1 = axes;
ax1.Visible = 'off';
set(gca,'visible','off')
imagesc(backgroundImagec);
colormap(ax1,'gray');
axis off;
ax2 = axes;
imagesc(ax2,heatmapImagec,'alphadata',heatmapTransparency);
colormap(ax2,jet); % or default 'jet' for higher contrast
%  delta1=(0/100)*(max(max(heatmapImage))-min(min(heatmapImage)));
%  delta2=(0/100)*(max(max(heatmapImage))-min(min(heatmapImage)));
%caxis(ax2,[min(min(heatmapImage))+delta1 max(max(heatmapImage))-delta2]);
caxis(ax2,[0 1]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
set(gca,'visible','off')
axis off;
colorbar;

hold on; plot(x3-50,y3-50,'color','r','linewidth',3);
drawnow; % Force it to draw immediately.
% plot(x2,y2, 'y','linewidth',3);
% drawnow; % Force it to draw immediately.
% plot(x,y, 'color','c','linewidth',3);
% drawnow; % Force it to draw immediately.
end