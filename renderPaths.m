function []=renderPaths(idxs)

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

folder = 'Data/Paths/';
pathfiles = dir(sprintf('%spath_*.dat', folder));
if nargin ~= 1 || idxs == -1
  idxs = 1:length(pathfiles);
  idxs = 1:128;
end
for f = idxs
    % cla
    if nargin == 1
      data = readBezier(sprintf('%spath_%d.dat', folder, f));
      fprintf('File: %s\n', sprintf('%spath_%d.dat', folder, f))
    else
      data = readBezier(sprintf('%s%s', folder, pathfiles(f).name));
      fprintf('File: %s\n', pathfiles(f).name)
    end
    
    for i = 1:data.n-1
        X = B0 * data.p(i,1) + B1 * data.out(i,1) + B2 * data.in(i+1,1) + B3 * data.p(i+1,1);
        Y = B0 * data.p(i,2) + B1 * data.out(i,2) + B2 * data.in(i+1,2) + B3 * data.p(i+1,2);
        plot(X, Y, 'color', [0.5, 0.5, 0.5]);
    end

    % drawnow
    % pause(.1)
    % fprintf('Path %d.\n', f);
end

fprintf('Files: %d\n', length(pathfiles))

plot(data.p(1,1),data.p(1,2),'bs','MarkerSize',8,'MarkerFaceColor','y');
plot(data.p(end,1),data.p(end,2),'kp','MarkerSize',10,'MarkerFaceColor','g');

axis equal
box on
grid on
xlim([-50 50])
ylim([-50 50])
xticks(-50:25:50)
yticks(-50:25:50)
