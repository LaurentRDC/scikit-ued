function [C,numberOfOverlapMaskedPixels] = normxcorr2_masked(varargin)

% [C,numberOfOverlapMaskedPixels] =
% normxcorr2_masked(fixedImage, movingImage, fixedMask, movingMask)
%   Masked normalized two-dimensional cross-correlation.
%   This function calculates the Masked NCC using FFTs
%   instead of spatial correlation.  It is therefore much faster for
%   larger structuring elements.
%   The masked normalized cross-correlation of fixedImage and
%   movingImage are calculated using masks fixedMask and movingMask.
%   fixedMask and movingMask should consist of only 1s and 0s, where 1
%   indicates locations of useful information in the corresponding image,
%   and 0 indicates locations that should be masked (ignored).
%   fixedImage and movingImage need not be the same size, but fixedMask
%   must be the same size as fixedImage, and movingMask must be the same
%   size as movingImage.
%   The resulting matrix C contains correlation coefficients and its values
%   range from -1.0 to 1.0. 
%
%   References: 
%   D. Padfield. "Masked Object Registration in the Fourier Domain".
%   Transactions on Image Processing. 
%   D. Padfield. "Masked FFT registration". In Proc. Computer Vision and
%   Pattern Recognition, 2010. 
%
%   Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
%

[fixedImage, movingImage, fixedMask, movingMask] = ParseInputs(varargin{:});
clear varargin;

fixedMask = double(fixedMask);
movingMask = double(movingMask);

% Ensure that the masks consist of only 0s and 1s.  Anything less than or
% equal to 0 gets set to 0, and everything else gets set to 1.
fixedMask(fixedMask <= 0 ) = 0;
fixedMask(fixedMask > 0) = 1;
movingMask(movingMask <= 0) = 0;
movingMask(movingMask > 0) = 1;

% The fixed and moving images need to be masked for the equations below to
% work correctly.
fixedImage(fixedMask==0) = 0;
movingImage(movingMask==0) = 0;

% Flip the moving image and mask in both dimensions so that its correlation
% can be more easily handled.
rotatedMovingImage = rot90(movingImage,2);
rotatedMovingMask = rot90(movingMask,2);
clear movingImage movingMask;

% Calculate all of the FFTs that will be needed.
fixedImageSize = size(fixedImage);
movingImageSize = size(rotatedMovingImage);
combinedSize = fixedImageSize + movingImageSize - 1;
% Find the next largest size that is a multiple of a combination of 2, 3,
% and/or 5.  This makes the FFT calculation much faster.
optimalSize(1) = FindClosestValidDimension(combinedSize(1));
optimalSize(2) = FindClosestValidDimension(combinedSize(2));

% Only 6 FFTs are needed.
fixedFFT = fft2(fixedImage,optimalSize(1),optimalSize(2));
rotatedMovingFFT = fft2(rotatedMovingImage,optimalSize(1),optimalSize(2));
fixedMaskFFT = fft2(fixedMask,optimalSize(1),optimalSize(2));
rotatedMovingMaskFFT = fft2(rotatedMovingMask,optimalSize(1),optimalSize(2));

% Only 6 IFFTs are needed.
% Compute and save these results rather than computing them multiple times.
numberOfOverlapMaskedPixels = real(ifft2(rotatedMovingMaskFFT.*fixedMaskFFT));
numberOfOverlapMaskedPixels = round(numberOfOverlapMaskedPixels);
numberOfOverlapMaskedPixels = max(numberOfOverlapMaskedPixels,eps);
maskCorrelatedFixedFFT = real(ifft2(rotatedMovingMaskFFT.*fixedFFT));
maskCorrelatedRotatedMovingFFT = real(ifft2(fixedMaskFFT.*rotatedMovingFFT));

numerator = real(ifft2(rotatedMovingFFT.*fixedFFT)) - ...
    maskCorrelatedFixedFFT .* maskCorrelatedRotatedMovingFFT ./ numberOfOverlapMaskedPixels;
clear rotatedMovingFFT fixedFFT

fixedSquaredFFT = fft2(fixedImage.*fixedImage,optimalSize(1),optimalSize(2));
fixedDenom = real(ifft2(rotatedMovingMaskFFT.*fixedSquaredFFT)) - ...
    maskCorrelatedFixedFFT.^2 ./ numberOfOverlapMaskedPixels;
fixedDenom = max(fixedDenom,0);
clear rotatedMovingMaskFFT fixedSquaredFFT maskCorrelatedFixedFFT;

rotatedMovingSquaredFFT = fft2(rotatedMovingImage.*rotatedMovingImage,optimalSize(1),optimalSize(2));
movingDenom = real(ifft2(fixedMaskFFT.*rotatedMovingSquaredFFT)) - ...
    maskCorrelatedRotatedMovingFFT.^2 ./ numberOfOverlapMaskedPixels;
movingDenom = max(movingDenom,0);
clear fixedMaskFFT rotatedMovingSquaredFFT maskCorrelatedRotatedMovingFFT;    

denom = sqrt(fixedDenom .* movingDenom);
clear fixedDenom movingDenom;

% denom is the sqrt of the product of positive numbers so it must be
% positive or zero.  Therefore, the only danger in dividing the numerator
% by the denominator is when dividing by zero.  
% Since the correlation value must be between -1 and 1, we therefore
% saturate at these values.
C = zeros(size(numerator));
tol = 1000*eps( max(abs(denom(:))) );
i_nonzero = find(denom > tol);
C(i_nonzero) = numerator(i_nonzero) ./ denom(i_nonzero);
C(C < -1) = -1;
C(C > 1) = 1;

% Crop out the correct size.
C = C(1:combinedSize(1),1:combinedSize(2));
numberOfOverlapMaskedPixels = numberOfOverlapMaskedPixels(1:combinedSize(1),1:combinedSize(2));

%-----------------------------------------------------------------------------
function [fixedImage, movingImage, fixedMask, movingMask] = ParseInputs(varargin)

iptchecknargin(4,4,nargin,mfilename)

fixedImage = varargin{1};
movingImage = varargin{2};
fixedMask = varargin{3};
movingMask = varargin{4};

iptcheckinput(fixedImage,{'logical','numeric'},{'real','nonsparse','2d','finite'},mfilename,'fixedImage',1)
iptcheckinput(movingImage,{'logical','numeric'},{'real','nonsparse','2d','finite'},mfilename,'movingImage',2)
iptcheckinput(fixedMask,{'logical','numeric'},{'real','nonsparse','2d','finite'},mfilename,'fixedMask',3)
iptcheckinput(movingMask,{'logical','numeric'},{'real','nonsparse','2d','finite'},mfilename,'movingMask',4)

% If either fixedImage or movingImage has a minimum value which is negative, we
% need to shift the array so all values are positive to ensure numerically
% robust results for the normalized cross-correlation.
fixedImage = shiftData(fixedImage);
movingImage = shiftData(movingImage);

%-----------------------------------------------------------------------------
function B = shiftData(A)

B = double(A);

is_unsigned = isa(A,'uint8') || isa(A,'uint16') || isa(A,'uint32');
if ~is_unsigned
    
    min_B = min(B(:)); 
    
    if min_B < 0
        B = B - min_B;
    end
    
end

%-----------------------------------------------------------------------------
function [newNumber] = FindClosestValidDimension(n)

% Find the closest valid dimension above the desired dimension.  This
% will be a combination of 2s, 3s, and 5s.

% Incrementally add 1 to the size until
% we reach a size that can be properly factored.
newNumber = n;
result = 0;
newNumber = newNumber - 1;
while( result ~= 1 )
    newNumber = newNumber + 1;
    result = FactorizeNumber(newNumber);
end

%-----------------------------------------------------------------------------
function [n] = FactorizeNumber(n)

for ifac = [2 3 5]
    while( rem(n,ifac) == 0 )
        n = n/ifac;
    end
end