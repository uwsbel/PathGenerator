function [] = renderBezier(data, np)
% Render the Bezier curve encoded in the 'data' structure using the
% specified number of points
%
% Radu Serban, 9/28/2107%

t = linspace(0, 1, np);
omt = 1 - t;

% Bernstein coefficients
B0 = omt.^3;
B1 = 3 .* t .* omt.^2;
B2 = 3 .* t.^2 .* omt;
B3 = t.^3;

% render curve segment between knots i and i+1
hold on
for i = 1 : data.n
    plot3(data.p(i,1), data.p(i,2), data.p(i,3), 'go');
end
for i = 1 : data.n-1
    X = B0 * data.p(i,1) + B1 * data.out(i,1) + B2 * data.in(i+1,1) + B3 * data.p(i+1,1);
    Y = B0 * data.p(i,2) + B1 * data.out(i,2) + B2 * data.in(i+1,2) + B3 * data.p(i+1,2);
    Z = B0 * data.p(i,3) + B1 * data.out(i,3) + B2 * data.in(i+1,3) + B3 * data.p(i+1,3);
    plot3(X, Y, Z, 'g-');
    
    % Uncomment to also plot control polygon 
    %plot3([data.p(i,1) data.out(i,1)], [data.p(i,2) data.out(i,2)], [data.p(i,3) data.out(i,3)], 'b-');
    %plot3([data.in(i+1,1) data.p(i+1,1)], [data.in(i+1,2) data.p(i+1,2)], [data.in(i+1,3) data.p(i+1,3)], 'r-');
end

