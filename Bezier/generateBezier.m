function data = generateBezier(p)
% Wrapper to generate a Bezier curve through given waypoints
%
% Radu Serban, 9/28/2107

data.n = size(p,1);
data.p = p;
[data.in, data.out] = splineBezier(p);
end
