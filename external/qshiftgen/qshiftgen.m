function [h0,g0,h1,g1,h,Ha,Dha,f0,T] = qshiftgen(param)

% function [h0,g0,h1,g1,h,Ha,Dha,f0,T] = qshiftgen(param)
% Copyright N G Kingsbury, Cambridge University, December 2002.
%
% Quarter-sample shift CWT design based on a linear-phase
% 2x oversampled even-length filter (length 2N).
% This version minimises the energy over a defined frequency range.
% It works with the minimum number of unknowns and specifies
% 2 predefined zeros at z=-1.
% This version uses dh2 as the update variable.
% This version uses fmin = 1/3 (= 1/6 of sampling rate) and minimises
% the energy of product of image filter gain and designed filter gain
% above f = 1/3. Above  f = 1/2, the image filter gain is assumed to be unity. 
% This version uses the same notation as the ICIP03 paper.
%
% Copyright Nick Kingsbury, University of Cambridge, January 2003.

% Test scripts:
% figure;N=10;h0 = qshiftgen([N 1/3 1 1 1]);title(sprintf('N = %d',N))
% figure;N=N+2;h0 = qshiftgen([N 1/3 1 1 1]);title(sprintf('N = %d',N))


if nargin < 1, param = []; end
if length(param) < 1, param(1) = 14; end
if length(param) < 2, param(2) = 1/3; end
if length(param) < 3, param(3) = 2; end
if length(param) < 4, param(4) = 1; end
if length(param) < 5, param(5) = 1; end

N = param(1); % No of filter taps (must be even, typ 12 to 18).
fmin = param(2); % Bottom of stopband (typ 1/3).
nzpi = param(3); % No of zeros at pi in the final filter (ie at pi/2 in h).
nzhpi = param(4); % No of zeros at pi in h (must be odd).
loopgain = param(5); % Gain for updating h in loop.
if rem(N,2)>0, 
   error('Wavegen: Qshift filter order must be even!');
end
if rem(nzhpi,2)~=1, 
   error('Wavegen: No of zeros at pi in h must be odd!');
end

% Design initial guess at prototype h filter.
np = 64;
mp = 4; 
x = [0:np]' * (1/np);
theta = [x.^mp; 2-x(np:-1:2).^mp] * (pi/4);
Hp = [cos(theta); zeros(2*np,1)];  % Prototype freq response.

% Apply half-sample phase shift and do IFFT.
Hp = Hp .* exp([0:(4*np-1)]' * (j*pi/(8*np)));
hp = real(ifft([Hp; 0; conj(Hp((4*np):-1:2))]));
% Truncate impulse reponse to 2N taps and ensure it is symmetric.
h = 2 * hp([N:-1:1  1:N]);

% Form zeros at pi in final filter, ie at pi/2 in h.
hpi = 1;
for k = 1:nzpi,
   hpi = conv(hpi,[1;0;1]);
end
% Form extra zeros at pi in h.
for k = 2:nzhpi,
   hpi = conv(hpi,[1;1]);
end

% Convolution matrix for hpi.
H = toeplitz([hpi; zeros(2*N-length(hpi),1)],[hpi(1) zeros(1,2*N-length(hpi))]); 
% Convert H to use only one side of symmetric error vector Dh.
H = H(:,1:end/2) + H(:,end:-1:end/2+1);

[h1,err_engy,Berr] = deconv_sym(h,hpi); %Deconvolve hpi from h.
h2 = h1(1:(end/2)); % h2 = first half of h1.
h = H*h2;
initgain=sum(h)
%  H0 = fft(h,8*np);

% Test PR accuracy of prototype.
%  pf=conv(h,h).'*100;
%  mat2scr([0 pf(4:4:4*(N-1))],'%9.4f');

% Loop to find h which satisfies PR
% and minimises energy of the frequency components defined by f0.
mag = 1;
f0 = pi*(fmin + (1-fmin)*[0:0.25:N].'/N);
if1 = find(f0 < 0.5*pi); % Find indices to frequencies less than pi/2.
% Fourier matrix for quadrupled frequency components to synthesise image of low passband.
f4 = exp(j*4*f0(if1)*([1:(length(h)/2)] - 1));
ss = 1:2:length(h);
T = ones(size(f0));
F = exp(j*f0*([1:length(h)] - 1));
c = [zeros(N/2-1,1); 1];
Ha = [];
Dha = [];
for i=1:40,
   Cfull = toeplitz([h; zeros(2*N-1,1)],[h(1) zeros(1,2*N-1)]); % Full convolution matrix.
   C = Cfull(4:4:2*N,:);  % Pick out every 4th row in C as far as the middle row.
	T(if1) = abs(f4 * (h(ss)/sum(h(ss)))); % Calculate freq response of image filter.
   TF = (T * ones(1,length(h))) .* F; % Scale F by image filter response, T.
   beta = min(2^((i-1)/1),2^20); % Increasing scale factor for PR terms.
   CFH = [2*beta*C; TF] * H;  % Complete LHS matrix.
   rhs = [beta*(c-C*h); -TF*h]; % RHS of equation.
   Dh = real(CFH\rhs); % Find Dh for best fit to PR and freq responses.
   h  = h + H*Dh*loopgain;  % Update h.
   Ha = [Ha abs(fft(h,512))];
   Dha = [Dha Dh];
end
finalgain=sum(h)
% figure(gcf);semilogy(rms(Dha)) % Plot convergence of Dh.

err = CFH*Dh - rhs;
rmserr1 = rms(err(1:N/2));
rmserr2 = rms(err(N/2+1:end));
h0 = h(2:2:2*N);
% figure;semilogy([0:256]*2/512,max(Ha(1:257,end),1e-4),'Linewidth',2);grid on
% title(sprintf('N = %d, nzpi = %d, nzhpi = %d',N,nzpi,nzhpi))
% hold on; semilogy(f0/pi,T,'r','Linewidth',2); hold off

p0 = conv(h0,rev(h0));
rmserrp = rms(p0(2:2:(end-1)/2));

h0pr = prforce(h0);  % Force PR on final result.
prerr = rms(h0pr-h0);
%h0 = h0pr;
fprintf('rmserr1 = %f, rmserr2 = %f, rmserrp = %g, prerr = %f\n',rmserr1,rmserr2,rmserrp,prerr);
h0 = sqrt(2) * h0 / sum(h0);
h0 = h0 / sum(h0);

g0 = h0(N:-1:1);   % g0 is time reverse of h0.
h1 = -g0 .* cumprod(-ones(size(g0)));  % h1 and g1 are quadrature mirrors
g1 = h0 .* cumprod(-ones(size(h0))); % of g0 and h0 respectively.

return;