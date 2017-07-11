% MaskedTranslationRegistrationTest
%   Test the MaskedTranslationRegistration code.
%
%   References: 
%   D. Padfield. "Masked Object Registration in the Fourier Domain".
%   Transactions on Image Processing. 
%   D. Padfield. "Masked FFT registration". In Proc. Computer Vision and
%   Pattern Recognition, 2010. 
%
%   Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
%


close all;
clear variables;
clc;

% Perform the registration on several sets of images.
x = [75 -130 130];
y = [75 130 130];
overlapRatio = 1/10;

for i = 1:3
    fixedImage = imread(sprintf('OriginalX%2iY%2i.png',x(i),y(i)));
    movingImage = imread(sprintf('TransformedX%2iY%2i.png',x(i),y(i)));
    fixedMask = fixedImage~=0;
    movingMask = movingImage~=0;

    [translation,maxC] = MaskedTranslationRegistration(fixedImage,movingImage,fixedMask,movingMask,overlapRatio);
    [transformedMovingImage] = TransformImageTranslation(movingImage,translation);

    % Given the transform, transform the moving image.
    [overlayImage] = OverlayRegistration(fixedImage,transformedMovingImage);
    figure; imagesc(overlayImage); title(['Test ' num2str(i) ': Registered Overlay Image']);

    disp(['Test ' num2str(i) ':']);
    disp(['Computed translation: ' num2str([translation(1) -translation(2)])]);
    disp(['Correlation score: ' num2str(maxC)]);
    trueTranslation = [x(i),-y(i)];
    disp(['Transformation error: ' num2str(translation - trueTranslation)]);
    disp(' ');
end
