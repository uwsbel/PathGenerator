function [in, out] = splineBezier(p)
% Generate Bezier curve through given waypoints
%
% Radu Serban, 9/28/2107

np = size(p,1);
in = zeros(np,3);
out = zeros(np,3);

in(1,:) = p(1,:);
out(np,:) = p(np,:);

% Special case for two points
if (np == 2)
   out(1,:) = (2*p(1,:) + p(2,:))/3;
   in(2,:) = (p(1,:) + 2*p(2,:))/3;
   return
end

n = np - 1;

% x coordinate
rhs = 4 * p(1:end-1,1) + 2 * p(2:end,1);
rhs(1) = p(1,1) + 2 * p(2,1);
rhs(end) = (8 * p(end-1,1) + p(end,1))/2;
x = solveTridiag(rhs, n);

% y coordinate
rhs = 4 * p(1:end-1,2) + 2 * p(2:end,2);
rhs(1) = p(1,2) + 2 * p(2,2);
rhs(end) = (8 * p(end-1,2) + p(end,2))/2;
y = solveTridiag(rhs, n);

% z coordinate
rhs = 4 * p(1:end-1,3) + 2 * p(2:end,3);
rhs(1) = p(1,3) + 2 * p(2,3);
rhs(end) = (8 * p(end-1,3) + p(end,3))/2;
z = solveTridiag(rhs, n);

for i = 0:n-2
    out(i+1,1) = x(i+1); 
    out(i+1,2) = y(i+1); 
    out(i+1,3) = z(i+1);
    
    in(i+2,1) = 2 * p(i+2,1) - x(i+2);
    in(i+2,2) = 2 * p(i+2,2) - y(i+2);
    in(i+2,3) = 2 * p(i+2,3) - z(i+2);
end
out(n,1) = x(n);
out(n,2) = y(n);
out(n,3) = z(n);

in(n+1,1) = (p(n+1,1) + x(n))/2;
in(n+1,2) = (p(n+1,2) + y(n))/2;
in(n+1,3) = (p(n+1,3) + z(n))/2;
end

function x = solveTridiag(rhs, n)
x = zeros(n,1);
tmp = zeros(n,1);

b = 2;
x(1) = rhs(1)/b;

for i = 1:n-1
    tmp(i+1)=1/b;
    if (i==n-1)
        b = 3.5 - tmp(i+1);
    else
        b = 4 - tmp(i+1);
    end
    x(i+1) = (rhs(i+1)-x(i))/b;
end
for i = 1:n-1
    x(n-i) = x(n-i) - tmp(n-i+1)*x(n-i+1);
end
end