function [finalRed, finalGreen, finalBlue, address]=readCzi() 
%       Written by Vardges Tserunyan for Wickersham Lab at MIT in 2017.
%       This function uses data obtained from the MATLAB Bioformats 5.5.2
%   bfopen() function to return red- and green-channel images from a .czi
%   file. Each image contains data from only one color channel (R or G), so 
%   without proper adjustment, it would appear as greyscale. 
%
%       INPUTS   - none
%       OUTPUTS  - finalRed (image from the red channel as a uint16 matrix)
%                  finalGreen (image from the green channel as a uint16 matrix)
%                  address (full path to the file)
%       REQUIRES - MATLAB Bioformats 5.5.2. Other versions might not
%                  necessarily work.
%
    %% Upload the file and extract the images
    %   Choose a .czi file, and use bfopen() to get a cell array. The array
    % contains matrices for each Z-plane of red or green color, as well as
    % metadata.
    addpath('/users/Vardges/Documents/MATLAB/bfmatlab_552')
    [filename, filepath]=uigetfile('*.czi');
    address=fullfile(filepath, filename);
    imageStack=bfopen(address);
    
    %% Separate channels
    %   Use the colormap to find out the color channel of a z-plane and add 
    % it to a correponding stack. Repeat for each of the planes
    for layer=1:length(imageStack{1,1})
    
        colorScheme=imageStack{1,3}{1,layer}(end,:);
        color=find(colorScheme);
        redStack(:,:,layer)=(imageStack{1,1}{layer,1})*uint16(color==1);
        greenStack(:,:,layer)=(imageStack{1,1}{layer,1})*uint16(color==2);
        blueStack(:,:,layer)=(imageStack{1,1}{layer,1})*uint16(color==3);
    end
    
    %% Create a Z-projection and adjust for brightness
    %   Having obtained a stack for both channels, take a Z-max projection
    % to cocentrate the stacks for each color channel into one image.
    finalRed=max(redStack,[],3);
    finalGreen=max(greenStack,[],3);
    finalBlue=max(blueStack,[],3);
    
    %   Adjust the brightness of the image to match a certain target
    % histogram. The histogram has been found empirically to provide the
    % most optimal brightness and is consistent for all images used. 
    load('parameters.mat', 'targetHgram');
    finalRed=histeq(finalRed,targetHgram);
    finalGreen=histeq(finalGreen,targetHgram);
    finalBlue=histeq(finalBlue,targetHgram);
    
    %   To eliminate background noise arising as an artifact of brightness
    % adjustment, use a 6x6 averaging filter.
    avg=(1/36)*ones(6);
    finalRed=imfilter(finalRed, avg);
    finalGreen=imfilter(finalGreen, avg);
    finalBlue=imfilter(finalBlue, avg);
end