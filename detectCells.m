function [processedFourier processedNoise processedThresh croppingMask]=detectCells(channel, color)
%       Written by Vardges Tserunyan for Wickersham Lab at MIT in 2017.
%       This function detects cells on an input imageby using predetermined 
% parameters. Green and redchannels are treated slightly differently: for a 
% green image, theuser gets to crop out the relevant fragement for analysis, 
% leaving out conspicuous sourses of error. For red images, nothing is 
% cropped: the function operates on the whole image. This difference comes 
% from the way overlap between green and red cells is analyzed at a later 
% stage of processing.
%
%       INPUTS  - channel (a uint16 image)
%                 color (can be only 1 for red or 2 for green)
%       OUTPUTS - processedFourier (the image after Fourier filtering)
%                 processedNoise (the image after morphological opening)
%                 processedThresh (the image after thresholding)
%                 croppingMask (the analyzed part of the image)
%

    %% Initial adjustments
    %   Load parameters that are used later for cell detection and crop out the
    % relevant section to analyze on a green image.
    load('parameters.mat', 'fourierRadius', 'diskRadius', 'thresh');
    croppingMask=uint16(ones(size(channel)));
    if color==2
        croppingMask=uint16(roipoly(channel));
    end
    close gcf
    %% Fourier filtering
    %   Pass the image through a high-pass Fourier filter to get rid of the
    % background making small-scale features (including cells) more distinct.
    highpass=zeros(size(channel));
    midpoint=round(size(channel)/2);
    [y,x]=size(channel);
    [y,x]=meshgrid(1:x,1:y);
    highpass=(sqrt((x-midpoint(1)).^2+(y-midpoint(2)).^2))>fourierRadius;

    processedFourier=fft2(channel);
    fourierShifted=fftshift(processedFourier);
    fourierShifted=fourierShifted.*highpass;
    processedFourier=(ifft2(ifftshift(fourierShifted)));
    processedFourier=uint16(real(processedFourier));

    %% Morphological opening
    %   Use morphological opening to eliminate small-scale noisy features, such
    % as cell processes.
    window=strel('disk', diskRadius, 0);
    processedNoise= imopen(processedFourier, window); 
    
    %% Segmentation
    %   Use a threshold to segment the cells.
    processedThresh= (processedNoise>=thresh*256).*double(croppingMask);
end
