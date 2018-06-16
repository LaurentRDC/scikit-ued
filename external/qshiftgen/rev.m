function y = rev(x)
% function y = rev(x)
% Reverse the order of the rows of x (ie invert each column).

y = x(size(x,1):-1:1,:);
return

