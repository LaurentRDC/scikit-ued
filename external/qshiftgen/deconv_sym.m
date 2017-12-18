function [y,err_engy,Berr] = deconv_sym(B,A)

% function [y,err_engy,Berr] = deconv_sym(B,A)
% Symmetrical deconvolution of A into B, such that B approx = conv(A,y)
% using matrix pseudo-inverse method.
% Berr is the error between B and conv(A,y).
% err_engy is the enrgy of Berr.
%
% Nick Kingsbury, Cambridge University, Dec 2001.

% Form a convolution matrix of A.
z = zeros(length(B) - length(A),1);
Cm = toeplitz([A(1); z],[A(:); z]).';

% Solve for y which gives LMS error for Cm * y = B
y = Cm \ B;

Berr = B - Cm * y;
err_engy = real(Berr' * Berr); 
return;