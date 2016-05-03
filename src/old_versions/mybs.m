function x = mybs(a,z)
%
%MYBS demo backward substitution:
% solves U*x = z
[m,n] = size(a);
x = zeros(n,1);
nm = min(n,m);
x(1:nm) = z(1:nm);
for i=nm:-1:1
    x(i) = x(i)/a(i,i);
    x(1:i-1) = x(1:i-1) - a(1:i-1,i)*x(i);
end