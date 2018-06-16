% QSHTEST2.M
% qshift wavelet freq response test M-file.
% Plots waveforms and frequency responses for complex wavelets at levels 1 to 6.
% param is the argument passed to qshiftgen.m .
%
% Use these lines of code to test this routine.
% param = [18 1/3 1 1 1]; qshtest2
% param = [24 1/3 1 1 1]; qshtest2
%
% Copyright Nick Kingsbury, University of Cambridge, January 2003.

% figure(1);
[ha,ga,h1a,g1a] = qshiftgen(param);
% sc = 1/sqrt(2);
% ha = ha * sc; ga = ga * sc; h1a = h1a * sc; g1a = g1a * sc;

% Use near-symmetric (13,19)-tap biorth. filters for level 1 wavelet filters.
% Could have used Daub (9,7) filters instead.
load near_sym_b
h0n = [h0o; 0];
g0n = [0; h0o];
g1n = [h1o; 0];
h1n = [0; h1o];

figure(5); 
cplot(beside(h0n+j*g0n,g1n+j*h1n),1),grid, pause
figure(6)
plot([-512:511]/256,max(abs(fftshift([fft(h0n+j*g0n,1024) ...
      fft(g1n+j*h1n,1024)])),1e-4))
axis([-2 2 0 2]),grid,pause
nf = 2048;
for i = 1:6,
   h1n = conv(h0n,upm(h1a,2^i));
   g1n = conv(g0n,upm(g1a,2^i));
   h0n = conv(h0n,upm(ha,2^i));
   g0n = conv(g0n,upm(ga,2^i));
   figure(5); cplot([h0n+j*g0n  g1n+j*h1n],2^(-i)),grid, pause
   figure(6+i)
%   plot([0:1023]/256,20*log10(max(abs([fft(h0n,1024)  fft(h1n,1024)]),1e-4)))
%   axis([0 4 -80 0]),grid,pause
plot([-nf:nf-1]*2/nf,max(abs(fftshift([fft(h0n+j*g0n,2*nf) ...
      fft(g1n+j*h1n,2*nf)])),1e-4))
   axis([-2 2 0 2]),grid,pause

end

