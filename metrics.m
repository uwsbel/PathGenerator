function [avg_dev, max_dev, path_length] = metrics(lfile, ffile, plotit)
% Calculate metrics for follower path deviation.
% Usage:
%   [avg_dev, max_dev, path_length] = metrics(lfile, ffile, plotit)
% Inputs:
%   lfile - filename (csv) with leader's path data
%   ffile - filename (csv) with follower's path data
%   plotit - if true, plot intermediate results
% Input files are assumed to have x-y coordinates in 2nd and 3rd column.
% Outputs:
%   avg_dev - average integrated deviation
%   max_dev - maximum deviation
%   path_length - length of leader path (only portion used in calculations)

addpath(genpath('SampleCurve/'));

% Load data for leader and follower
ldata=csvread(lfile,1,0);
fdata=csvread(ffile,1,0);

% Leader and follower paths (x-y)
% (Eliminate duplicates else spline complains)
lpath = unique(ldata(:,2:3),'stable','rows');
fpath = unique(fdata(:,2:3),'stable','rows');

% Assuming the follower is always behind the leader, calculate:
% - p_start = closest point on follower path to the leadder's start point
% - p_end   = closest point on leader path to the follower's end point
p_start = distance2curve(fpath, lpath(1,:), 'spline');
p_end = distance2curve(lpath, fpath(end,:), 'spline');

if plotit
    figure('position', [200 200 800 800])
    hold on
    hp_l = plot(lpath(:,1),lpath(:,2));
    hp_f = plot(fpath(:,1),fpath(:,2));
    hp_s = plot(p_start(1),p_start(2),'o', 'markerfacecolor', 'g');
    hp_e = plot(p_end(1),p_end(2),'o', 'markerfacecolor', 'r');
end

% Clamp paths to a "common" range.
% - Find the closest point on follower path to p_start. Discard all points before that
% - Find closest point on leader path to p_end. Discard all points after that
[~,ii] = min( sum((fpath - p_start).^2, 2) );
fpath = fpath(ii:end,:);
[~,ii] = min( sum((lpath - p_end).^2, 2) );
lpath = lpath(1:ii,:);

% Estimate length of leader's path
%%path_length = sum( vecnorm(diff(lpath),2,2) );
%%fprintf('leader path length = %f\n', path_length);

% Evaluate the leader path at 'n' equally spaced points in arclength
n = min(500, size(lpath,1));
lpath = interparc(n,lpath(:,1),lpath(:,2),'spline');

if plotit
    hp_l2 = plot(lpath(:,1),lpath(:,2),'.');
    axis equal
    box on
    grid on
    xlim([-50 50])
    ylim([-50 50])
    xticks(-50:25:50)
    yticks(-50:25:50)
    
    legend([hp_l hp_f hp_s hp_e hp_l2], ... 
           { ...
             'Leader path (original)', ...
             'Follower path (original)', ...
             'Clip point (start)', ...
             'Clip point (end)',...
             'Leader path (clipped & resampled)' ...
           }, ...
           'location', 'southeast')
end

% Estimate arclength between two consecutive points and total length of
% leader's path
ii = round(n/2);
delta = vecnorm(lpath(ii,:)-lpath(ii-1,:), 2, 2);
path_length = (n-1) * delta;
%%fprintf('leader path length = %f\n', path_length);

% For each point on the leader's path, find closest distance to follower
% path and accumulate metrics.  The deviation is seen as a function of
% cummulative distance along the leader's path.
tic
[~,deviation,~] = distance2curve(fpath, lpath, 'spline');
avg_dev = sum(deviation) / n;
max_dev = max(deviation);
toc

if plotit
   figure
   plot((1:n)*delta, deviation)
   grid on, box on
   xlabel('Distance along leader path [m]')
   ylabel('Follower path deviation [m]')
end
