function y=beside(x1,x2)

% function y=beside(x1,x2)
% Put x1 beside x2, padding with zeros.

[m1,n1]=size(x1);
[m2,n2]=size(x2);
m = max(m1,m2);
y=zeros(m,n1+n2);
y((m-m1)/2+[1:m1],[1:n1])=x1;
y((m-m2)/2+[1:m2],n1+[1:n2])=x2;
return
