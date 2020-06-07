function [xx,yy] = evaluateBezier(XS, YS, t)
% Generate a Bezier curve through given waypoints and evaluate at
% specified points
%
% Radu Serban, 9/28/2107

data = generateBezier([XS' YS' zeros(length(XS),1)]);

omt = 1 - t;

% Bernstein coefficients
B0 = omt.^3;
B1 = 3 .* t .* omt.^2;
B2 = 3 .* t.^2 .* omt;
B3 = t.^3;

% Evaluate points on Bezier curve
for i = 1 : data.n-1
    xx = B0 * data.p(i,1) + B1 * data.out(i,1) + B2 * data.in(i+1,1) + B3 * data.p(i+1,1);
    yy = B0 * data.p(i,2) + B1 * data.out(i,2) + B2 * data.in(i+1,2) + B3 * data.p(i+1,2);
end
