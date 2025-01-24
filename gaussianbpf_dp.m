function filtered_image = gaussianbpf_dp(I,d0,d1)
%% Function adjusted by DPouli from original file exchange function from :
% Leonardo O. Iheme (leonardo.iheme@cc.emu.edu.tr)
% 24th of March 2011
%
%   I = The input grey scale image
%   d0 = Lower cut off frequency
%   d1 = Higher cut off frequency
%
% The function makes use of the simple principle that a bandpass filter
% can be obtained by multiplying a lowpass filter with a highpass filter
% where the lowpass filter has a higher cut off frquency than the high pass filter.
%
% Usage GAUSSIANBPF(I,DO,D1)
% Example
% ima = imread('grass.jpg');
% ima = rgb2gray(ima);
% filtered_image = gaussianbpf(ima,30,120);
% Gaussian Bandpass Filter


%%

f = double(I);
[nx ny] = size(f);
% f = uint8(f);
fftI = fft2(f,2*nx-1,2*ny-1);
fftI = fftshift(fftI);

% figure;
% subplot(1,2,1)
% imshow(f,[]);
% title('Original Image')

%subplot(2,2,2)
%fftshow(fftI,'log')
%title('Fourier Spectrum of Image')
% Initialize filter.
filter1 = ones(2*nx-1,2*ny-1);
filter2 = ones(2*nx-1,2*ny-1);
filter3 = ones(2*nx-1,2*ny-1);

for i = 1:2*nx-1
    for j =1:2*ny-1
        dist = ((i-(nx+1))^2 + (j-(ny+1))^2)^.5;
        % Use Gaussian filter.
        filter1(i,j) = exp(-dist^2/(2*d1^2)); % higher cutt off (cuts everything above that)
        filter2(i,j) = exp(-dist^2/(2*d0^2)); % lower cut off ( cuts everything above that)
        filter3(i,j) = 1.0 - filter2(i,j); % invert lower cut off so it cuts everything below that
        filter3(i,j) = filter1(i,j).*filter3(i,j); %multiply 1st and 3rd to get your bandpass result
    end
end
% Update image with passed frequencies
filtered_image = fftI + filter3.*fftI;

%subplot(2,2,3)
%fftshow(filter3,'log')
%title('Frequency Domain Filter Function Image')
filtered_image = ifftshift(filtered_image);
filtered_image = ifft2(filtered_image,2*nx-1,2*ny-1);
filtered_image = real(filtered_image(1:nx,1:ny)); % recover edited image
% filtered_image = uint8(filtered_image);

% subplot(1,2,2)
% imshow(filtered_image,[])
%title('Bandpass Filtered Image')