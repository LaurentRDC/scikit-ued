function [transformedImage] = TransformImageTranslation(movingImage,transform,fixedImageSize)

% [transformedImage] = TransformImageTranslation(movingImage,transform,fixedImageSize)
%   Transform the moving image to line up with the fixed image to enable
%   visual comparison of the registered images. 
%   movingImage is the image to transform using the two-element vector
%   transform.  If fixedImageSize is not specified, it will be set to the size of
%   movingImage.
%   Note that a positive transform will shift the image to the left and up,
%   whereas a negative transform will shift the image to the right and down.
%   This is because the transform is applied to the indices that are then
%   used to sample the image pixels.  This is consistent with standard
%   transform equations.
%   Intuitively, this is like shifting the bounding box of the image
%   according to the transform and then setting the image to be the pixels
%   under the transformed bounding box.
%   This is also consistent with the ITK definitions.
%
%   References: 
%   D. Padfield. "Masked Object Registration in the Fourier Domain".
%   Transactions on Image Processing. 
%   D. Padfield. "Masked FFT registration". In Proc. Computer Vision and
%   Pattern Recognition, 2010. 
%
%   Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
%

if( nargin < 3 )
    fixedImageSize = size(movingImage);
end

movingImageSize = size(movingImage);

% Note that the brackets are necessary to keep the size of x the same as
% imageSize(2).
x = [1:fixedImageSize(2)] + transform(1);
borderIndicesX = find(x<1 | x > movingImageSize(2));
nonBorderIndicesX = find(x>=1 & x <= movingImageSize(2));
x(borderIndicesX) = 1;
y = [1:fixedImageSize(1)] + transform(2);
borderIndicesY = find(y<1 | y > movingImageSize(1));
nonBorderIndicesY = find(y>=1 & y <= movingImageSize(1));
y(borderIndicesY) = 1;

transformedImage = zeros(fixedImageSize);
transformedImage = movingImage(y,x);
transformedImage(:,borderIndicesX) = 0;
transformedImage(borderIndicesY,:) = 0;