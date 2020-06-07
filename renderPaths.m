addpath(genpath('Bezier/'));

t = linspace(0, 1, 100);
omt = 1 - t;
% Bernstein coefficients
B0 = omt.^3;
B1 = 3 .* t .* omt.^2;
B2 = 3 .* t.^2 .* omt;
B3 = t.^3;

figure
hold on

pathfiles = dir('path_*.dat');
for f = 1:length(pathfiles)
    data = readBezier(pathfiles(f).name);
    for i = 1:data.n-1
        X = B0 * data.p(i,1) + B1 * data.out(i,1) + B2 * data.in(i+1,1) + B3 * data.p(i+1,1);
        Y = B0 * data.p(i,2) + B1 * data.out(i,2) + B2 * data.in(i+1,2) + B3 * data.p(i+1,2);
        plot(X, Y, 'color', [0.5, 0.5, 0.5]);
    end
end

plot(data.p(1,1),data.p(1,2),'bs','MarkerSize',8,'MarkerFaceColor','y');
plot(data.p(end,1),data.p(end,2),'kp','MarkerSize',10,'MarkerFaceColor','g');

axis equal
box on
grid on
xlim([-50 50])
ylim([-50 50])
xticks(-50:25:50)
yticks(-50:25:50)
