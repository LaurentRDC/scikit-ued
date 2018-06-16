function y = upm(x,m)

% function y = upm(x,m)
% Up-sample the column vector x m:1 by inserting m-1 zeroes between each
% sample but exclude the final zeroes.
% If m is not given, m = 2.

if nargin < 2, m = 2; end

temp = [x zeros(size(x).*[1 m-1])].';
y = temp(:);
y = y(1:(length(y)-m+1));
return
