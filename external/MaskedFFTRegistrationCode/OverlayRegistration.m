function [outImage] = OverlayRegistration(fixedImage,transformedMovingImage)

% [outImage] = OverlayRegistration(fixedImage,transformedImage)
%   Overlay a fixed image with a transformed moving image.
%   The intensities of outImage are rescaled to uint8.
%
%   References: 
%   D. Padfield. "Masked Object Registration in the Fourier Domain".
%   Transactions on Image Processing. 
%   D. Padfield. "Masked FFT registration". In Proc. Computer Vision and
%   Pattern Recognition, 2010. 
%
%   Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
%

fixedImage = double(fixedImage);
transformedMovingImage = double(transformedMovingImage);

fixedImageSize = size(fixedImage);
movingImageSize = size(transformedMovingImage);

if( fixedImageSize(1) > movingImageSize(1) )
    borderWidth = [(fixedImageSize(1) - movingImageSize(1)) 0];
    transformedMovingImage = padarray(transformedMovingImage,borderWidth, 'post');
elseif( fixedImageSize(1) < movingImageSize(1) )
    borderWidth = [(movingImageSize(1) - fixedImageSize(1)) 0];
    fixedImage = padarray(fixedImage,borderWidth, 'post');
end

if( fixedImageSize(2) > movingImageSize(2) )
    borderWidth = [0 (fixedImageSize(2) - movingImageSize(2))];
    transformedMovingImage = padarray(transformedMovingImage,borderWidth, 'post');
elseif( fixedImageSize(2) < movingImageSize(2) )
    borderWidth = [0 (movingImageSize(2) - fixedImageSize(2))];
    fixedImage = padarray(fixedImage,borderWidth, 'post');
end

outImage = zeros(size(fixedImage,1),size(fixedImage,2),3,'uint8');

% Scale the images to 8 bit
fixedImage = double(fixedImage);
fixedImage = 255*(fixedImage-min(fixedImage(:)))/(max(fixedImage(:)) - min(fixedImage(:)));
transformedMovingImage = double(transformedMovingImage);
transformedMovingImage = 255*(transformedMovingImage-min(transformedMovingImage(:)))/(max(transformedMovingImage(:)) - min(transformedMovingImage(:)));

outImage(:,:,1) = uint8(fixedImage);
outImage(:,:,2) = uint8(transformedMovingImage);
outImage(:,:,3) = uint8(zeros(size(fixedImage)));