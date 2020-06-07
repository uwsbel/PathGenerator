function data = readBezier(file)
% Read Bezier curve from file in Chrono::Vehicle format
%
% Radu Serban, 9/28/2107

A = dlmread(file);
data.n = A(1,1);
nc = A(1,2);

if (nc == 9)
    data.p = A(2:end,1:3);
    data.in = A(2:end,4:6);
    data.out = A(2:end,7:9);
end

if (nc == 3)
    data.p = A(2:end,1:3);
    [data.in, data.out] = splineBezier(data.p);
end
