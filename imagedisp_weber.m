close all;
clc;
clear;
names={'1749','2616','2801','4390','4279','1450','4928','4978','3425','4155','4551','1341','1311','3936','3889','3109','2972','2891','2645','2347','2283'};

for i = 1:1%length(names)%patient2:1-35, patient4:36- ,patient6:

xi=51;xm=1000;yi=51;ym=1100;
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

rd = im2double(rd); % red diff
gd = im2double(gd);
bd = im2double(bd);
rs = im2double(rs); % red sum
gs = im2double(gs);
bs = im2double(bs);
rp=  im2double(rp); % red perpendicular
gp=  im2double(gp);
bp=  im2double(bp);

rsn=(rs-min(min(rs)))./(max(max(rs))-min(min(rs)));
gsn=(gs-min(min(gs)))./(max(max(gs))-min(min(gs)));
bsn=(bs-min(min(bs)))./(max(max(bs))-min(min(bs)));
rdn=(rd-min(min(rd)))./(max(max(rd))-min(min(rd)));
gdn=(gd-min(min(gd)))./(max(max(gd))-min(min(gd)));
bdn=(bd-min(min(bd)))./(max(max(bd))-min(min(bd)));
rpn=(rp-min(min(rp)))./(max(max(rp))-min(min(rp)));
gpn=(gp-min(min(gp)))./(max(max(gp))-min(min(gp)));
bpn=(bp-min(min(bp)))./(max(max(bp))-min(min(bp)));

% Set bilateral filter parameters.
w     = 5;       % bilateral filter half-width
sigma = [3 0.1]; % bilateral filter standard deviations

% Apply bilateral filter to each image.
% rdt = bfilter2(rdn,w,sigma); % red diff
% gdt = bfilter2(gdn,w,sigma);
% bdt = bfilter2(bdn,w,sigma);
% rst = bfilter2(rsn,w,sigma); % red sum
% gst = bfilter2(gsn,w,sigma);
% bst = bfilter2(bsn,w,sigma);
% rpt = bfilter2(rpn,w,sigma); % red perpendicular
% gpt = bfilter2(gpn,w,sigma);
% bpt = bfilter2(bpn,w,sigma); 

rd2 = glare(rdn); % red diff
gd2 = glare(gdn);
bd2 = glare(bdn);
rs2 = glare(rsn); % red sum
gs2 = glare(gsn);
bs2 = glare(bsn);
rp2 = glare(rpn); % red perpendicular
gp2 = glare(gpn);
bp2 = glare(bpn);

var11=im2double((rp2-gp2));
var12=(var11-min(min(var11)))./(max(max(var11))-min(min(var11)));
var12=(imcomplement(var12));
var1=(var12-min(min(var12)))./(max(max(var12))-min(min(var12)));
var21=im2double((rd2-gd2));
var22=(var21-min(min(var21)))./(max(max(var21))-min(min(var21)));
var22=(imcomplement(var22));
var2=(var22-min(min(var22)))./(max(max(var22))-min(min(var22)));

var3=im2double(imadjust((rp2./gp2)));
var41=im2double(rd2./gd2);
var4=im2double(imcomplement(var41));

% backgroundImage = imadjust(rd+gd+bd);
% backgroundImagec = backgroundImage(51:1070,51:1100);
% % percentage1 = 0;
% % percentage2 =1;
% % low_bound = percentage1 / 100;
% % high_bound = 1 - percentage2/100;
% % K = imadjust(Im1,stretchlim(Im1,[low_bound high_bound]),[]);
% 
% %heatmapImage1 = imgaussfilt(var1,20);
heatmapImage1 = gaussianbpf_dp(var11,25,256);
heatmapImage1 = gaussianbpf_dp(heatmapImage1,20,120);
heatmapImage1 = butterworthbpf_caroline(heatmapImage1,10,120,3);
heatmapImage1 = medfilt2(heatmapImage1);
% 
% %heatmapImage2 = imgaussfilt(var21,20);
% heatmapImage2 = gaussianbpf_dp(var21,25,256);
% heatmapImage2 = gaussianbpf_dp(heatmapImage2,20,120);
% heatmapImage2 = butterworthbpf_caroline(heatmapImage2,10,120,3);
% heatmapImage2 = medfilt2(heatmapImage2);
% 
% %heatmapImage3 = imgaussfilt(var3,20);
% heatmapImage3 = gaussianbpf_dp(var3,25,256);
% heatmapImage3 = gaussianbpf_dp(heatmapImage3,20,120);
% heatmapImage3 = butterworthbpf_caroline(heatmapImage3,10,120,3);
% heatmapImage3 = medfilt2(heatmapImage3);
% 
% %heatmapImage4 = imgaussfilt(var41,4);
% heatmapImage4 = gaussianbpf_dp(var41,25,256);
% heatmapImage4 = gaussianbpf_dp(heatmapImage4,20,120);
% heatmapImage4 = butterworthbpf_caroline(heatmapImage4,10,120,3);
% heatmapImage4 = medfilt2(heatmapImage4);
% 
% heatmapImage5=imadjust(im2double(0.2*heatmapImage1+0.4*heatmapImage2+0.8*heatmapImage3));
% %heatmapImage5=(im2double(rs+gs+bs));
% % mask=heatmapImage5(:,:)<=1;
% % heatmapImage=heatmapImage5.*mask;
% %  heatmapImage=(heatmapImage5);
% %heatmapImage=imadjust(im2double(rd+gd+bd));
% heatmapImage=(heatmapImage5-min(min(heatmapImage5)))./(max(max(heatmapImage5))-min(min(heatmapImage5)));
% heatmapImagec= heatmapImage5(51:1070,51:1100);

v1=heatmapImage1;
v2=heatmapImage1;
v1(~bwL1)=0;
v2(~bwL3)=0;
L1(i,1)=mean2(nonzeros(v1));L1(i,2)=std(nonzeros(v1));
L3(i,1)=mean2(nonzeros(v2));L3(i,2)=std(nonzeros(v2));
%weber(i)=(L1-L3)./L3
display(L1);display(L3);

% sizeHeatmapImage = size(heatmapImagec);
% sizeBackgroundImage = size(backgroundImagec);
% 
% heatmapTransparency = ones(size(backgroundImagec,1),size(backgroundImagec,2));
% heatmapTransparency = heatmapTransparency*0.6;
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

% %images sum and diff
% bw=bwL3(xi:xm,yi:ym);
% bound = bwboundaries(bw,'noholes');
% n=size(bound,1);
% hold on;
% for k=1:n
%     thisbound=bound{k};
%     x=thisbound(:,2);
%     y=thisbound(:,1);
%     plot(x,y,'color','r','linewidth',3);
%     drawnow;
% end
% figure;rgb=cat(3,rst,gst,bst);rgb=rgb(xi:xm,yi:ym,:);imshow(rgb);
% hold on;
% for k=1:n
%     thisbound=bound{k};
%     x=thisbound(:,2);
%     y=thisbound(:,1);
%     plot(x,y,'color','r','linewidth',3);
%     drawnow;
% end
% figure;rgb=cat(3,rdt,gdt,bdt);rgb=rgb(xi:xm,yi:ym,:);imshow(rgb);
% hold on;
% for k=1:n
%     thisbound=bound{k};
%     x=thisbound(:,2);
%     y=thisbound(:,1);
%     plot(x,y,'color','r','linewidth',3);
%     drawnow;
% end
% figure;rgb=cat(3,rpt,gpt,bpt);rgb=rgb(xi:xm,yi:ym,:);imshow(rgb);
% hold on;
% for k=1:n
%     thisbound=bound{k};
%     x=thisbound(:,2);
%     y=thisbound(:,1);
%     plot(x,y,'color','r','linewidth',3);
%     drawnow;
% end
end