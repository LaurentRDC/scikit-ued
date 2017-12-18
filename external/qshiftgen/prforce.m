function hpr = prforce(h)

% function hpr = prforce(h)
% Solve for symmetric ladder filter to implement the filter h
% and its time-reverse quadrature mirror h1.
% Return the value of the perfect reconstruction version of h.
% h must be a column vector.
%
% Nick Kingsbury, Cambridge University, Feb 99.

n = length(h);
n2 = fix(n/2) - 1;

% Form the highpass filter to go with h, by reversing h and inverting
% the sign of alternate terms.
h1 = h(n:-1:1) .* cumprod(-ones(n,1));
hh=[h h1];
a = []; 
d = [];
err = [];

% Loop for each stage in the ladder to find the ladder coefficent
% a(i) which most accurately cancels out the first two elements of
% hh(:,1) and the last two elements of hh(:,2) or vice-versa,
% depending on which gives the lower energy result.
for i = 1:n2,
	sh = size(hh);
	t = 3:sh(1);
	a(i) = sum(prod(hh(1:2,:).')) / sum(hh(1:2,2).^2); % Find least squares solution for a(i).
   if abs(a(i)) < 1,
      hh1 = hh * [1 a(i); -a(i) 1]; % Do inverse ladder filtering.
      hh = [hh1(t,1)  hh1(t-2,2)]; % Trim off the terms in hh1 which should be approx zero.
      d(i) = 1; % d(i) defines which decision was taken.
   else
      a(i) = 1/a(i);
      hh1 = hh * [1 -a(i); a(i) 1];
      hh = [hh1(t-2,1)  hh1(t,2)];
      d(i) = 2;
   end
   err(:,i) = hh1(1:2,d(i)); % Save the approximation errors.
end

% hh2 = hh
% a,d,err

% Reconstruct hh from the a(i) values, using the d(i) also.
z2 = [0;0];
for i = n2:-1:1,
	sh = size(hh);
	t = 3:sh(1);
   if d(i) == 1,
      hh1 = [[z2;hh(:,1)]  [hh(:,2);z2]];
      hh = hh1 * [1 -a(i); a(i) 1]; % Do forward ladder filtering
   else
      hh1 = [[hh(:,1);z2]  [z2;hh(:,2)]];
      hh = hh1* [1 a(i); -a(i) 1];
   end
end

% Correct for the gain change due to inverse and forward ladder filtering.
hpr = hh(:,1) / prod(1+a.*a);
return

