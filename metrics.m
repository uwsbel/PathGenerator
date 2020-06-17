function [mPath, mSpeed, mSpeedRel, path_length] = metrics(lfile, ffile, varargin)
% Calculate metrics for follower path deviation.
% Usage:
%   [mPath, mSpeed, mSpeedRel, path_length] = metrics(lfile, ffile)
%   [mPath, mSpeed, mSpeedRel, path_length] = metrics(lfile, ffile, plotit)
% Inputs:
%   lfile - filename (csv) with leader's path data
%   ffile - filename (csv) with follower's path data
%   timeit - if true, report timing information
%   plotit - if true, plot intermediate results
% Input files are assumed to have x-y coordinates in 2nd and 3rd column.
% Outputs:
%   mPath - metrics for path deviation
%     mPath.avg - average integrated deviation
%     mPath.max - maximum deviation
%   mSpeed - metrics for speed deviation
%     mSpeed.avg - average integrated deviation
%     mSpeed.max - maximum deviation
%   path_length - length of leader path (only portion used in calculations)

if nargin > 2
    plotit = varargin{1};
else
    plotit = false;
end

addpath(genpath('SampleCurve/'));

tstart_total = tic;

% Load data for leader and follower
ldata=csvread(lfile,1,0);
fdata=csvread(ffile,1,0);

% Leader and follower paths (x-y)
% (Eliminate duplicates else spline complains)
[lpath, lidx] = unique(ldata(:,2:3),'stable','rows');
[fpath, fidx] = unique(fdata(:,2:3),'stable','rows');

% Leader and follower horizontal speed
lspeed = vecnorm(ldata(lidx, 5:6), 2, 2);
fspeed = vecnorm(fdata(fidx, 5:6), 2, 2);

clear lidx fidx

% Assuming the follower is always behind the leader, calculate:
% - p_start = closest point on follower path to the leadder's start point
% - p_end   = closest point on leader path to the follower's end point
p_start = distance2curve(fpath, lpath(1,:), 'spline');
p_end = distance2curve(lpath, fpath(end,:), 'spline');

% Clip paths to a "common" range.
% - Find the closest point on follower path to p_start. Discard all points before that
% - Find closest point on leader path to p_end. Discard all points after that
[~,fidx] = min( sum((fpath - p_start).^2, 2) );
fpathC = fpath(fidx:end,:);
[~,lidx] = min( sum((lpath - p_end).^2, 2) );
lpathC = lpath(1:lidx,:);

% Clip speeds to the same interval
fspeedC = fspeed(fidx:end);
lspeedC = lspeed(1:lidx);

clear fidx lidx

% Estimate length of leader's path
%%path_length = sum( vecnorm(diff(lpath),2,2) );
%%fprintf('leader path length = %f\n', path_length);

% Evaluate the leader path at 'n' equally spaced points in arclength
n = min(200, size(lpathC,1));
tstart_sampling = tic;
lpathS = interparc(n, lpathC(:,1), lpathC(:,2), 'spline');
time_sampling = toc(tstart_sampling);

% Estimate arclength between two consecutive points and total length of
% leader's path
ii = round(n/2);
delta = vecnorm(lpathS(ii,:)-lpathS(ii-1,:), 2, 2);
path_length = (n-1) * delta;
%%fprintf('leader path length = %f\n', path_length);

% For each point on the leader's path, find closest point on follower path
% and corresponding distances (these are the follower's deviations).
tstart_closest = tic;
[fpathS,deviation,~] = distance2curve(fpathC, lpathS, 'spline');
time_closest = toc(tstart_closest);

% Calculate deviation metrics.  The deviation is seen as a function of
% cummulative distance along the leader's path.
mPath.avg = sum(deviation) / n;
mPath.max = max(deviation);

% Find vehicle speeds at each sampled point on the leader path and at the
% corresponding closest point on the follower path. We approximate by
% using the speed from the closest point on a clipped path (leader or follower)
% to each sample point.
tstart_speed = tic;
lspeedS = zeros(n,1);
fspeedS = zeros(n,1);
for k = 1:n
[~, lidx] = min(vecnorm(lpathC - lpathS(k,:), 2, 2));
lspeedS(k) = lspeedC(lidx);
[~, fidx] = min(vecnorm(fpathC - fpathS(k,:), 2, 2));
fspeedS(k) = fspeedC(fidx);
end
time_speed = toc(tstart_speed);

% Note: for speed metric calculations, discard an initial segment of the
% common path (always a fixed value).  This is because:
% - experiments atsrt with all vehicles at rest
% - they all need to accelerate from rest to reach the "convoy speed"
% - at t=0, the relative speed deviation would involve a division by 0
% Ideally, the test setup would be modified to collect data between to
% pre-defined "start" and "end" points on the prescribed path with all
% followers using PID controllers before reaching the "start" point and
% then switching to an inference controller.
buf_length = 10;  % meters
is = round(buf_length / delta);
ns = n - is + 1;
fspeedS = fspeedS(is:end);
lspeedS = lspeedS(is:end);


% Calculate speed metrics (absolute speed deviation)
speed_deviation = abs(fspeedS - lspeedS);
mSpeed.avg = sum(speed_deviation) / ns;
mSpeed.max = max(speed_deviation);

% Calculate speed metrics (absolute speed deviation)
speed_deviation_rel = speed_deviation ./ abs(lspeedS);
mSpeedRel.avg = sum(speed_deviation_rel) / ns;
mSpeedRel.max = max(speed_deviation_rel);

clear lidx fidx

time_total = toc(tstart_total);

fprintf('Time: %.4f  %.4f  %.4f  %.4f\n', time_sampling, time_closest, time_speed, time_total);

if plotit
    % Colors for leader and follower
    lcol = [0, 0.4470, 0.7410];
    fcol = [0.85, 0.325, 0.098];
    
    % Figure 1:  path calculations
    figure('position', [200 200 800 800])
    hold on
    
    % Input leader and follower paths
    hp_l = plot(lpath(:,1),lpath(:,2), 'color', lcol);
    hp_f = plot(fpath(:,1),fpath(:,2), 'color', fcol);
    
    % Clip points (start on follower path, end on leader path)
    hp_s = plot(p_start(1),p_start(2),'o', 'color', fcol, 'markerfacecolor', fcol);
    hp_e = plot(p_end(1),p_end(2),'o', 'color', lcol, 'markerfacecolor', lcol);
    
    % Sampled leader path (equal arclength)
    hp_l2 = plot(lpathS(:,1),lpathS(:,2),'.', 'color', lcol);

    % Closest points on follower path
    hp_f2 = plot(fpathS(:,1),fpathS(:,2),'.', 'color', fcol);
    
    % Deviations
    for i = 1:n
       plot([lpathS(i,1) fpathS(i,1)], [lpathS(i,2), fpathS(i,2)], 'color', [0.5, 0.5, 0.5], 'linewidth', 0.5) 
    end

    axis equal
    box on
    grid on
    xlim([-50 50])
    ylim([-50 50])
    xticks(-50:25:50)
    yticks(-50:25:50)
    title('Vehicle paths')

    legend([hp_l hp_f hp_s hp_e hp_l2 hp_f2], ...
        { ...
        'Leader path (original)', ...
        'Follower path (original)', ...
        'Clip point (start)', ...
        'Clip point (end)',...
        'Leader path (clipped & resampled)' ...
        'Follower path (closest points)' ...
        }, ...
        'location', 'southeast')
    
    % Figure 2: Deviation as function of distance traveled along leader path
    figure
    plot((1:n)*delta, deviation)
    grid on, box on
    xlabel('Distance along leader path [m]')
    ylabel('Follower path deviation [m]')
    title('Path deviation')
    
    % Figure 3: Vehicle speeds at same locations along leader path
    figure
    plot((is:n)*delta, lspeedS, 'color', lcol)
    hold on, grid on, box on
    plot((is:n)*delta, fspeedS, 'color', fcol)
    xlabel('Distance along leader path [m]')
    ylabel('Vehicle speed [m/s]')
    title('Vehicle speeds at same locations')
    legend('Leader', 'Follower', 'location', 'southeast')
    
    % Figure 4: Speed deviation as function of distance traveled along leader path
    figure
    plot((is:n)*delta, speed_deviation)
    grid on, box on
    xlabel('Distance along leader path [m]')
    ylabel('Follower speed deviation [m/s]')
    title('Speed deviation')

    % Figure 5: Speed deviation as function of distance traveled along leader path
    figure
    plot((is:n)*delta, speed_deviation_rel)
    grid on, box on
    xlabel('Distance along leader path [m]')
    ylabel('Follower speed deviation [-]')
    title('Relative speed deviation')

end
