%   Written by Vardges Tserunyan for Wickersham lab at MIT in 2017.
%   This script loads a given .czi image through redCzi() and then detects
% cells through detectCells(). Then, it goes on to count individual green
% cells first and, secondly, only the ones that also appear on the red
% image. The script returns these numbers, as well as a .tif  file with
% images detailing the process of analysis.

%   ****
%   INSTALL BIOFORMATS 5.5.2 BEFORE USE
%   ****

%% Finding the cells and returning images
%   Load a .czi file and detect cells on green and red images
[channelRed, greenChannel, ~,address]=readCzi();
[filepath, filename, ~]=fileparts(address);
[greenFourier, greenNoise, greenThresh, croppingMask]=detectCells(greenChannel, 2);
[redFourier, redNoise, ~, ~]=detectCells(channelRed, 1);

%% Detecting, counting and returning overlap fraction
%   For each green cell, detect the region that includes the cell body and 
% another region that both includes and surrounds the cell. Then, go to the 
% red image and see if the correposnding region in it contains a cell 
% or not. Assume a region contains a red cell, if its mean intensity is 
% significantly higher than the mean of the region including the surround.
cellLocsGreen= bwconncomp(greenThresh);
cellSurrGreen= bwconncomp(bwmorph(greenThresh, 'thicken', 5));

avgCells=regionprops(cellLocsGreen, redNoise, 'MeanIntensity');
avgEnv=regionprops(cellSurrGreen, redNoise, 'MeanIntensity');
intensityDiff=cell2mat(struct2cell(avgCells))-cell2mat(struct2cell(avgEnv));
cellLocsGreen.Overlaps=intensityDiff>50;
numOverlap=sum(cellLocsGreen.Overlaps);

fprintf('Number of green cells: %f\n', cellLocsGreen.NumObjects)
fprintf('Number of overlaps: %f\n', numOverlap)
fprintf('The fraction of overlap is %f\n', numOverlap/cellLocsGreen.NumObjects)
filename
%% Separating green-and-red cells from green-not-red cells
%   Create a matrix that contains the outlines of all green cells that also
% appear in the red channel.
positiveIndices=find(cellLocsGreen.Overlaps);
positiveMatrix=zeros(size(channelRed));
for k=1:length(positiveIndices)
    positiveMatrix=positiveMatrix+(labelmatrix(cellLocsGreen)==positiveIndices(k));
end
outlinedPositive=bwperim(positiveMatrix);
outlinedPositive=bwmorph(outlinedPositive, 'thicken', 1);
outlinedPositive=im2uint16(outlinedPositive);

%   Create a matrix that contains the outlines of all green cells that do
% not appear in the red channel.
negativeMatrix=(labelmatrix(cellLocsGreen)>0) - double(positiveMatrix);
outlinedNegative=bwperim(negativeMatrix);
outlinedNegative=bwmorph(outlinedNegative, 'thicken', 1);
outlinedNegative=im2uint16(outlinedNegative);

%% Saving processed images to a .tif file
%   Create a .tif file and save it in the same folder under the original 
% image name plus "_detected.tif" in the end. The .tif file contains, in the 
% following order (1)the brightness-adjusted green channel image from the .czi 
% file, (2) alldetected cells on the green image, outlined as red circles if  
% they alsoappear on the red channel and blue otherwise, (3) the brightness-adjusted 
% red channel image from the .czi file, (4) Green cells that do not appear 
% in the red image outlined, (5) Green cells that do appear on the red
% image outlined.
outlinedNegativeIm=cat(3, zeros(size(greenChannel)), greenChannel.*croppingMask, outlinedNegative);
outlinedPositiveIm=cat(3, outlinedPositive, greenChannel.*croppingMask, zeros(size(greenChannel)));
outlinedGreenAll=cat(3, outlinedPositive, greenChannel.*croppingMask, outlinedNegative);
greenChannel=cat(3,zeros(size(greenChannel)), greenChannel, zeros(size(greenChannel)));
channelRed=cat(3,channelRed, zeros(size(channelRed)), zeros(size(channelRed)));

imwrite(greenChannel, [address '_' num2str(cellLocsGreen.NumObjects) '_' num2str(numOverlap) '.tif'  ])
imwrite(outlinedGreenAll, [address '_' num2str(cellLocsGreen.NumObjects) '_' num2str(numOverlap) '.tif' ], 'WriteMode', 'append')
imwrite(channelRed, [address '_' num2str(cellLocsGreen.NumObjects) '_' num2str(numOverlap)  '.tif'], 'WriteMode', 'append')
imwrite(outlinedNegativeIm, [address '_' num2str(cellLocsGreen.NumObjects) '_' num2str(numOverlap)  '.tif'], 'WriteMode', 'append')
imwrite(outlinedPositiveIm, [address '_' num2str(cellLocsGreen.NumObjects) '_' num2str(numOverlap) '.tif'], 'WriteMode', 'append')
imwrite(double(croppingMask), [address '_' num2str(cellLocsGreen.NumObjects) '_' num2str(numOverlap) '.tif'], 'WriteMode', 'append')

% imwrite(greenChannel, 'OriginalGreen.tif')
% imwrite(cat(3, zeros(size(greenNoise)), greenFourier, zeros(size(greenNoise))), 'GreenFourier.tif')
% imwrite(cat(3, zeros(size(greenNoise)), greenNoise, zeros(size(greenNoise))), 'GreenOpening.tif')
% imwrite(double(croppingMask), 'CroppingMask.tif')
% imwrite(channelRed, 'OriginalRed.tif')
% imwrite(cat(3, zeros(size(greenNoise)), greenThresh, zeros(size(greenNoise))), 'GreenThresh.tif')
% imwrite(cat(3, redFourier, zeros(size(greenNoise)), zeros(size(greenNoise))), 'RedFourier.tif')
% imwrite(cat(3, redNoise, zeros(size(greenNoise)), zeros(size(greenNoise))), 'RedOpening.tif')
% imwrite(outlinedGreenAll, 'Final.tif')
% imwrite(cat(3, channelRed(:,:,1), greenChannel(:,:,2), zeros(size(greenChannel(:,:,1)))), 'Yellow.tif')

% subplot(5,2,1)
% imshow(greenChannel)
% title('Original green image', 'FontSize', 5)
% 
% subplot(5,2,2)
% imshow(channelRed)
% title('Original red image', 'FontSize', 5)
% 
% subplot(5,2,3)
% imshow(cat(3, zeros(size(greenNoise)), greenFourier, zeros(size(greenNoise))))
% title('Green image after Fourier high-pass filter', 'FontSize', 5)
% 
% subplot(5,2,4)
% imshow(cat(3, redFourier, zeros(size(greenNoise)), zeros(size(greenNoise))))
% title('Red image after Fourier high-pass filter', 'FontSize', 5)
% 
% subplot(5,2,5)
% imshow(cat(3, zeros(size(greenNoise)), greenNoise, zeros(size(greenNoise))))
% title('Green image after morphological opening', 'FontSize', 5)
% 
% subplot(5,2,6)
% imshow(cat(3, redNoise, zeros(size(greenNoise)), zeros(size(greenNoise))))
% title('Red image after morphogical opening', 'FontSize', 5)
% 
% subplot(5,2,7)
% imshow(double(croppingMask))
% title('Cropping mask to analyze the relavant ROI', 'FontSize', 5)
% 
% subplot(5,2,8)
% imshow(cat(3, zeros(size(greenNoise)), greenThresh, zeros(size(greenNoise))))
% title('Green image after thresholding and selecting the ROI', 'FontSize', 5)
% 
% subplot(5,2,9)
% imshow(outlinedGreenAll)
% title('Final result after overlap detection', 'FontSize', 5)
% 
% subplot(5,2,10)
% imshow(cat(3, channelRed(:,:,1), greenChannel(:,:,2), zeros(size(greenChannel(:,:,1)))))
% title('Two originl images overlaid', 'FontSize', 5)

