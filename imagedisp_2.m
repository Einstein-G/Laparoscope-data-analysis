close all;
clc;
clear;
warning('off');
names={'1749','2616','2751','4390','4279','1450','4928','4978','4978','3425','4155','4551','1341','1311','3936','3889','3109','2972','2891','2645','2347','2283','0023','0037','2615','2769','3321','3388','2265','3631','0064','0503','0432','1413','1848'};% 2801 replaced with 2751
names2={'1749','2616','2751','4390','4279','1450','4928','4978a','4978b','3425','4155','4551','1341','1311','3936','3889','3109','2972','2891','2645','2347','2283','0023','0037','2615','2769','3321','3388','2265','3631','0064','0503','0432','1413','1848'};% 2801 replaced with 2751
j=0;
fprintf('\t\tLesion\t\t\t\tHeatlthy\n');
for i = 30:30%length(names)%patient2:1-35, patient4:36- ,patient6:

%xi=51;xm=1050;yi=51;ym=1150;
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
rd = (im2double(rd)); % red diff
gd = (im2double(gd));
bd = (im2double(bd));
rs = (im2double(rs)); % red sum
gs = (im2double(gs));
bs = (im2double(bs));
rp=  (im2double(rp)); % red perpendicular
gp=  (im2double(gp));
bp=  (im2double(bp));

diff = cat(3,rd,gd,bd);
sum = cat(3,rs,gs,bs);
perp = cat(3,rp, gp, bp);

diff=im2double(diff);
diff_lab = rgb2lab(diff);
max_luminosity = 100;
L = diff_lab(:,:,1)./max_luminosity;
%L=glare(L);
L = imadjust(L);
%L = wiener2(adapthisteq(L),[5 5]);
%L = (adapthisteq(L));
L=imgaussfilt(L,2);

% L = gaussianbpf_dp(L,25,256);
% L = gaussianbpf_dp(L,20,120);
% L = butterworthbpf_caroline(L,10,120,3);
L=((L-min(min(L(xi:xm,yi:ym))))./(max(max(L(xi:xm,yi:ym)))-min(min(L(xi:xm,yi:ym)))))*max_luminosity;
diff_lab(:,:,1)=L;
diff1 = lab2rgb(diff_lab);
figure; imshow(diff1(xi:xm,yi:ym,:));
%figure; imshow(diff);
end
% 
% rd2 = glare(rd); % red diff
% gd2 = glare(gd);
% bd2 = glare(bd);
% rs2 = glare(rs); % red sum
% gs2 = glare(gs);
% bs2 = glare(bs);
% rp2 = glare(rp); % red perpendicular
% gp2 = glare(gp);
% bp2 = glare(bp);
% 
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
% 
% var11=normalise(double((w1*rp2+w2*gp2+w3*bp2+w4*rd2+w5*gd2+w6*bd2)));
% %vart(~bwL4)=0;
% backgroundImage = normalise(imadjust(rd+gd+bd));
% backgroundImagec = backgroundImage(xi:xm,yi:ym);
% 
% 
% heatmapImage1 = imgaussfilt(var11,20);
% 
% heatmapImage=normalise(imadjust(heatmapImage1));
% g=imgaussfilt(heatmapImage,30).*g;
% heatmapImage(g~=0)=g(g~=0);
% heatmapImagec= heatmapImage(xi:xm,yi:ym);
% 
% sizeHeatmapImage = size(heatmapImagec);
% sizeBackgroundImage = size(backgroundImagec);
% 
% heatmapTransparency = ones(size(backgroundImagec,1),size(backgroundImagec,2));
% heatmapTransparency = heatmapTransparency*1;
% %heatmapTransparency(~bwL4)=0;
% 
% % Showing background image and overlaying heatmap
% figure;
% ax1 = axes;
% ax1.Visible = 'off';
% set(gca,'visible','off')
% imagesc(backgroundImagec);
% colormap(ax1,'gray');
% axis off;
% ax2 = axes;
% imagesc(ax2,heatmapImagec,'alphadata',heatmapTransparency);
% colormap(ax2,jet); % or default 'jet' for higher contrast
% %  delta1=(0/100)*(max(max(heatmapImage))-min(min(heatmapImage)));
% %  delta2=(0/100)*(max(max(heatmapImage))-min(min(heatmapImage)));
% %caxis(ax2,[min(min(heatmapImage))+delta1 max(max(heatmapImage))-delta2]);
% caxis(ax2,[0 1]);
% ax2.Visible = 'off';
% linkprop([ax1 ax2],'Position');
% set(gca,'visible','off')
% axis off;
% colorbar;
% 
% bw1=bwL1(xi:xm,yi:ym);
% bound = bwboundaries(bw1,'noholes');
% n=size(bound,1);
% hold on;
% for k=1:n
%     thisbound=bound{k};
%     x=thisbound(:,2);
%     y=thisbound(:,1);
%     polyshapeL1 = polyshape(x,y);
% areaL1 = area(polyshapeL1);
% [x1,y1] = boundary(polyshapeL1);
% [cx1, cy1] = centroid(polyshapeL1); %find the centroid
% 
% n = 1;
% polyshapeL2 = 0;
% areaL2 = 0;
% while (areaL2 - areaL1) < (1.1*areaL1)
%     polyshapeL2 = polybuffer(polyshapeL1, n); %expand by n units
%     areaL2 = area(polyshapeL2);
%     n = n+0.1;
% end
% [x2,y2] = boundary(polyshapeL2);
% 
% % Calculating L3
% n = 1;
% polyshapeL3 = 0;
% areaL3 = 0;
% while (areaL3 - areaL2) < (1*areaL1)
%     polyshapeL3 = polybuffer(polyshapeL2, n); %expand by n units
%     areaL3 = area(polyshapeL3);
%     n = n+0.1;
% end
% [x3,y3] = boundary(polyshapeL3);
% 
% bwL1 = poly2mask(x1, y1, size(rs, 1), size(rs, 2)); %returns a binary mask size as input img
% bwL2 = poly2mask(x2, y2, size(rs, 1), size(rs, 2)); %returns a binary mask size as input img
% bwL3 = poly2mask(x3, y3, size(rs, 1), size(rs, 2)); %returns a binary mask size as input img
% bwL4 = bwL3;
% bwL3(bwL2)=0;
% 
%     plot(x3,y3,'color','r','linewidth',3);
%     drawnow;
%     temp1=heatmapImage;temp3=heatmapImage;
%     temp1(~bwL1)=0;
%     temp3(~bwL3)=0;
%     j=j+1;
%     L1(j,1)=mean2(nonzeros(temp1));L1(j,2)=std(nonzeros(temp1));
%     L3(j,1)=mean2(nonzeros(temp3));L3(j,2)=std(nonzeros(temp3));
%     L(j)=L1(j,1);
%     if(j<10)
%         gr(j)=1;
%     else
%         gr(j)=2;
%     end
%     fprintf('%f %s %f ',L1(j,1),char(177),L1(j,2));
%     fprintf('%f %s %f\n',L3(j,1),char(177),L3(j,2));
%     weber(j)=(L1(j,1)-L3(j,1))./L3(j,1);
% end
% end
% L(j+1:2*j)=L3(:,1);
% gr(j+1:2*j)=3;
% figure;boxplot(L,gr);
% 
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