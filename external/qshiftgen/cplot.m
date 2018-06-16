function cplot(h,sep)

% function cplot(h,sep)
% Plot a complex impulse response h in terms of its real and
% imag parts in red and blue.
% The magnitude is plotted in green if h is complex.
% Sep controls the vertical separation between the plot of each
% column of h.

if nargin < 2, sep = 0; end

[m,n]=size(h);
t = [1:m]' - (m+1)/2;

sepm = sep * ones(m,1) * [0:(n-1)];

plot(t, imag(h) - sepm, 'b-');
hold on
plot(t, real(h) - sepm, 'r-');
if ~isreal(h), 
  plot(t, abs(h) - sepm, 'g-');
end
hold off

return

