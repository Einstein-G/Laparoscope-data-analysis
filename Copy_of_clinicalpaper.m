close all;
clc;
clear;
warning('off');
names={'1749','2616','2751','4390','4279','1450','4928','4978','3425','4155','4551','1341','1311','3936','3889','3109','2972','2891','2645','2347','2283','0023','0037','2615','2769','3321','3388','2265','3631','0064','0503','0432','1413','1848'};% 2801 replaced with 2751
names2={'1749','2616','2751','4390','4279','1450','4928','4978','3425','4155','4551','1341','1311','3936','3889','3109','2972','2891','2645','2347','2283','0023','0037','2615','2769','3321','3388','2265','3631','0064','0503','0432','1413','1848'};% 2801 replaced with 2751
j=0;q=0;
B=[-13.2523441355976;22.0436432367643;-78.1882801588861;105.949466357600];
fprintf('\t\tLesion\t\t\t\tHeatlthy\n');
for i = 30:30%34:34%30:30%slength(names)%patient2:1-35, patient4:36- ,patient6:


filename=names{i};
matfilename=sprintf('%s%s%s','t0',names2{i},'.mat');
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

n1(i)=mean2(sum); n2(i)=mean2(diff);

%figure; imshow(sum);

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

w     = 5;       % bilateral filter half-width
%sigma = [3 0.1]; % bilateral filter standard deviations
sigma = [3 0.1]; % bilateral filter standard deviations

% rd2 = glare(rd); % red diff
% gd2 = glare(gd);
% bd2 = glare(bd);
% rs = glare(rs); % red sum
% gs = glare(gs);
% bs = glare(bs);
% rp = glare(rp); % red perpendicular
% gp = glare(gp);
% bp = glare(bp);

% rd = bfilter2(rd,w,sigma);
% gd = bfilter2(gd,w,sigma);
% bd = bfilter2(bd,w,sigma);
rs = bfilter2(rs,w,sigma);
gs = bfilter2(gs,w,sigma);
bs = bfilter2(bs,w,sigma);
% rp = bfilter2(rp,w,sigma);
% gp = bfilter2(gp,w,sigma);
% bp = bfilter2(bp,w,sigma);


diff = (cat(3,rd,gd,bd));
sum = cat(3,rs,gs,bs);
perp = cat(3,rp, gp, bp);


image=im2double(sum);
image_lab = rgb2lab(image);
max_luminosity = 100;
L = (image_lab(:,:,1)./max_luminosity);
LInv=imcomplement(L);
L = imcomplement(imreducehaze(LInv,'ContrastEnhancement','none'));
meanV = mean2(L);
% L = (meanV + 1.5 * (L - meanV));% Increase contrast by factor of 1.5
 L = L*6; L = min(L,1);
image_lab(:,:,1)   = (L).* 100;
 image_lab(:,:,2:3) = image_lab(:,:,2:3)*8; % Increase saturation
% L = (adapthisteq(L));
% L=glare(L);
%L = imadjust(L);
% L = bfilter2(L,w,sigma);
%L = wiener2(adapthisteq(L),[5 5]);

%L=imgaussfilt(L,2);

% L = gaussianbpf_dp(L,25,256);
% L = gaussianbpf_dp(L,20,120);
% L = butterworthbpf_caroline(L,10,120,3);
% L = (normalise(L)*max_luminosity);
% image_lab(:,:,1)=L;
%image_lab(:,:,2)=normalise(image_lab(:,:,2));
%image_lab(:,:,3)=normalise(image_lab(:,:,3));
image = imguidedfilter(lab2rgb(image_lab));

% hsvImage=rgb2hsv(image);
% hChannel = ((hsvImage(:, :, 1)));
% sChannel = (((hsvImage(:, :, 2)))); %sChannel(sChannel>1)=1;sChannel(sChannel<0)=0;
% vChannel = (normalise(bfilter2(imadjust(hsvImage(:, :, 3)),w,sigma)));
% 
% %vChannel = normalise(imgaussfilt((hsvImage(:, :, 3)),3));
% %meanV = mean2(vChannel);
% %newV = meanV + 1.5 * (vChannel - meanV); % Increase contrast by factor of 1.5
% newHSVImage = cat(3, hChannel, sChannel, vChannel);
% image = hsv2rgb(newHSVImage);

% rd = (im2double(rd)); % red diff
% gd = (im2double(gd));
% bd = (im2double(bd));
% rs = (im2double(rs)); % red sum
% gs = (im2double(gs));
% bs = (im2double(bs));
% rp=  (im2double(rp)); % red perpendicular
% gp=  (im2double(gp));
% bp=  (im2double(bp));
% rd2 = normalise(im2double(glare(rd))); % red diff
% gd2 = normalise(im2double(glare(gd)));
% bd2 = normalise(im2double(glare(bd)));
% rs2 = glare(rs); % red sum
% gs2 = glare(gs);
% bs2 = glare(bs);
% rp2 = glare(rp); % red perpendicular
% gp2 = glare(gp);
% bp2 = glare(bp);

% Set bilateral filter parameters.


% Apply bilateral filter to each image.
% rd3 = bfilter2(rd2,w,sigma); % red diff
% gd3 = bfilter2(gd2,w,sigma);
% bd3 = bfilter2(bd2,w,sigma);

% rd2 = imadjust(rd2,stretchlim(rd2,[0 0.99]));
% gd2 = imadjust(gd2,stretchlim(gd2,[0 0.99]));
% bd2 = imadjust(bd2,stretchlim(bd2,[0 0.99]));

%rd2 = gaussianbpf_dp(rd2,25,256);
%rd2 = gaussianbpf_dp(rd2,20,120);
%rd2 = butterworthbpf_caroline(rd,10,120,3);

%gd2 = gaussianbpf_dp(gd2,25,256);
%gd2 = gaussianbpf_dp(gd2,20,120);
%gd2 = butterworthbpf_caroline(gd,10,120,3);

%bd2 = gaussianbpf_dp(bd2,25,256);
%bd2 = gaussianbpf_dp(bd2,20,120);
%bd2 = butterworthbpf_caroline(bd,10,120,3);
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
w=1;
% % var11=normalise(im2double(((21.1155*rp2-5.6013*gp2-2.5904*bp2)-(bd2+gd2+rd2))));
% % %var11=normalise(im2double(((bp2+gp2+rp2)-(bd2+gd2+rd2))));
% % %var11=normalise(im2double(((bd2+gd2+rd2))));
% % %var11=normalise(im2double(rp2-(mean2(rp2)/mean2(bp2))*bp2));
% % var12=normalise(im2double((bd2+gd2+rd2)));
% % var13=normalise(im2double((var12-var11)));
% % var14=normalise(im2double(((mean(bd2(:))+mean(gd2(:))+mean(rd2(:)))/(mean(bp2(:))+mean(gp2(:))+mean(rp2(:))))*(bp2+gp2+rp2)-(bd2+gd2+rd2)));
% 
% 
% %results
% var21=normalise(double(((mean(rd2(:))/(mean(bd2(:))+mean(gd2(:))+mean(rd2(:))))*(bp2+gp2+rp2)-rd2)));
% var22=normalise(double(((mean(rd2(:))/(mean(bd2(:))+mean(gd2(:))+mean(rd2(:))))*(bp2+gp2+rp2)-(bd2+gd2+rd2))));
% var23=normalise(im2double(((bp2+gp2+rp2)-(bd2+gd2+rd2))));
% 
% % var1=normalise((imcomplement(var12)));
% % var2=normalise(imcomplement(var21));
% % var3=normalise(double(imadjust((rp2./bp2))));
% % var41=normalise(double(rd2./bd2));
% % var4=normalise(double(imcomplement(var41)));
% 
% %vart(~bwL4)=0;
% 
% backgroundImage = normalise(imadjust(rd+gd+bd));
% backgroundImagec = backgroundImage(xi:xm,yi:ym);
% 
% % w     = 50;       % bilateral filter half-width5
% % sigma = [50 0.5]; % bilateral filter standard deviations3 0.1
% 
% % Apply bilateral filter to each image.
% % rdt = bfilter2(rd,w,sigma); % red diff
% 
% heatmapImage1 = im2double(imgaussfilt(var11,3));
% % heatmapImage1 = var11;
% heatmapImage2 = im2double(imgaussfilt(var12,3));
% heatmapImage3 = im2double(heatmapImage2-heatmapImage1);
% %heatmapImage3 = imgaussfilt(var3,3);
% %heatmapImage4 = imgaussfilt(var41,3);
% %heatmapImage5=imadjust(im2double(0.2*heatmapImage1+0.4*heatmapImage2+0.8*heatmapImage3));
% 
% % change percentage to customize limits
% percentage1 = 2;%0
% percentage2 = 0;%0
% low_bound = percentage1 / 100;
% high_bound = 1 - percentage2/100;
% heatmapImage = imadjust(heatmapImage1,stretchlim(heatmapImage1,[low_bound high_bound]),[]);
% heatmapImage=normalise((im2double((heatmapImage))));
% 
% %heatmapImage=normalise((im2double(rs+gs+bs)));
% 
% % w     = 5;       % bilateral filter half-width
% % sigma = [3 0.1]; % bilateral filter standard deviations
% % 
% % % Apply bilateral filter to each image.
% % heatmapImage = bfilter2(heatmapImage,w,sigma); % red diff
% 
% % g=imgaussfilt(heatmapImage,30).*g;
% % heatmapImage(g~=0)=g(g~=0);
% heatmapImagec= heatmapImage(xi:xm,yi:ym);
% 
% sizeHeatmapImage = size(heatmapImagec);
% sizeBackgroundImage = size(backgroundImagec);
% 
%  heatmapTransparency = ones(size(backgroundImagec,1),size(backgroundImagec,2));
%  heatmapTransparency = heatmapTransparency*1;
% %heatmapTransparency(~bwL4)=0;
% % diffc=diff(xi:xm,yi:ym,:);
% % sumc=sum(xi:xm,yi:ym,:);
 imagec=image(xi:xm,yi:ym,:);
%imagec=image(1:400,300:1100,:);
% gsc=gd3(xi:xm,yi:ym);
% bsc=bd3(xi:xm,yi:ym);
% sumc = cat(3,rsc,gsc,bsc);
figure; imshow(imagec);
% Showing background image and overlaying heatmap
% figure;
% ax1 = axes;
% ax1.Visible = 'off';
% set(gca,'visible','off')
% imagesc(backgroundImagec);
% colormap(ax1,'gray');
% axis off;
clearvars sum; 
heatmapImage=(sum(image,3)./3);
heatmapImagec=heatmapImage(xi:xm,yi:ym,:);
figure;
ax2 = axes;
%imagesc(ax2,heatmapImagec,'alphadata',heatmapTransparency);
imagesc(ax2,heatmapImagec);
colormap(ax2,jet); % or default 'jet' for higher contrast
%  delta1=(0/100)*(max(max(heatmapImage))-min(min(heatmapImage)));
%  delta2=(0/100)*(max(max(heatmapImage))-min(min(heatmapImage)));
%caxis(ax2,[min(min(heatmapImage))+delta1 max(max(heatmapImage))-delta2]);
caxis(ax2,[0 1]);
ax2.Visible = 'off';
%linkprop([ax1 ax2],'Position');
set(gca,'visible','off')
axis off;
colorbar;

% bw1=logical(bwL1(xi:xm,yi:ym));
% bw3=logical(bwL3(xi:xm,yi:ym));
% temp3=heatmapImagec;
% temp3(bw3)=0;
%     q=q+1;
%     L3(q,1)=mean2(nonzeros(temp3));L3(q,2)=std(nonzeros(temp3));
% bound = bwboundaries(bw1,'noholes');
% n=size(bound,1);
% hold on;
% for k=2:n
%     thisbound=bound{k};
%     x=thisbound(:,2);
%     y=thisbound(:,1);
%     polyshapeL1 = polyshape(x,y);
% areaL1 = area(polyshapeL1);
%  [x1,y1] = boundary(polyshapeL1);
% % [cx1, cy1] = centroid(polyshapeL1); %find the centroid
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
% bwn1 = poly2mask(x1, y1, size(heatmapImagec, 1), size(heatmapImagec, 2)); %returns a binary mask size as input img
% bwn2 = poly2mask(x2, y2, size(heatmapImagec, 1), size(heatmapImagec, 2)); %returns a binary mask size as input img
% bwn3 = poly2mask(x3, y3, size(heatmapImagec, 1), size(heatmapImagec, 2)); %returns a binary mask size as input img
% bwn4 = bwn3;
% bwn3(bwn2)=0;
% 
% 
% %     plot(x3,y3,'color','r','linewidth',3);
% %     drawnow;
%     temp1=heatmapImagec;%temp3=heatmapImagec;
%     temp1(~bwn1)=0;
%     %temp3(bwn2)=0;
%     j=j+1;
%     L1(j,1)=mean2(nonzeros(temp1));L1(j,2)=std(nonzeros(temp1));
%     %L3(j,1)=mean2(nonzeros(temp3));L3(j,2)=std(nonzeros(temp3));
%     L(j)=L1(j,1);
%     if(j<10)
%         gr(j)=1;
%     else
%         gr(j)=2;
%     end
%     
% %     fprintf('%f %s %f ',L1(j,1),char(177),L1(j,2));
% 
% %    weber(j)=((L1(j,1)-L3(j,1))./L3(j,1));
% %    fprintf('%f\n',weber(j));
% end
end
% grn=gr';
% %weber=weber';
% %figure;boxplot(weber,grn,'Labels',{'malignant','Benign'});
% L(j+1:j+q)=L3(:,1);
% gr(j+1:j+q)=3;
% L=L';gr=gr';
% figure;boxplot(L,gr,'Labels',{'malignant','Benign','Healthy'});
% %ylim([0 1]);
% Lt= [L(1); L(2); L(4); L(56); L(65);L(66); L(67); L(69); L(121); L(130)];
% grt=[1; 1; 1; 2; 2;3; 3; 3; 3; 3];
% figure;boxplot(Lt,grt,'Labels',{'malignant','Benign','Healthy'});
% ylim([0 1]);

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