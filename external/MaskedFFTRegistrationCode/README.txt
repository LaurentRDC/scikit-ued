a. Description: These files include Matlab code for Masked FFT translation registration along with a test and several test images from the TIP paper "Masked Object Registration in the Fourier Domain" by Dirk Padfield.  The code can be used to register two images of interest that have associated masks.  If you use this code for your research, please cite this paper.

b. Size: 300kb

c. Platform: These files can be run using Matlab on any platform.

d. Environment: Matlab is required to run this code.

e. Major Component Description: Files for computing the translation component of the masked FFT registration along with tests and test images.

f. Detailed Set-up Instructions: Simply copy the files to a directory on your computer and point Matlab to it.

g. Detailed Run Instructions:
1) To run the test, simply type "MaskedTranslationRegistrationTest" on the Matlab command line.
2) To test masked FFT registration on your own images, pass them and their masks to the function "MaskedTranslationRegistration".
3) For more details, please see the help information for each function.

h. Output Description: The test will print the computed translation, the correlation score, and the transformation error for three pairs of images.  It will also display an overlay of each of the registered pairs.

i. Contact Information: Dirk Padfield, GE Global Research, padfield@research.ge.com