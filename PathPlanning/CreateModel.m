%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPAP115
% Project Title: Path Planning using PSO in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function model=CreateModel()
    % Domain
    xmin = -47;
    ymin = -47;
    xmax = 47;
    ymax = 47;

    % Start
    xs = -47;
    ys = -47;
    
    % Target (Destination)
    xt = 47;
    yt = 47;
    
    % ---------
    % Obstacles
    % ---------
    
    % Fixed (equally spaced between start and target)
    nfixed = 4;                         %  CAN BE CHANGED 
    foo = linspace(xs, xt, nfixed + 2);
    xfixed = foo(2:end-1);
    foo = linspace(ys, yt, nfixed + 2);
    yfixed = foo(2:end-1);
    rfixed = 6 * ones(1, nfixed);       %  CAN BE CHANGED
    
    % Random
    nrnd = 40;                          %  CAN BE CHANGED
    
    rmin = 2;                           %  CAN BE CHANGED
    rmax = 6;                           %  CAN BE CHANGED
    
    xrnd = (xmax - xmin - 8) .* rand(1,nrnd) + (xmin + 4);
    yrnd = (ymax - ymin - 8) .* rand(1,nrnd) + (ymin + 4);
    rrnd = (rmax - rmin) .* rand(1,nrnd) + rmin;
    
    % All obstacles
    xobs = [xfixed xrnd];
    yobs = [yfixed yrnd];
    robs = [rfixed rrnd];

    % ===================================================================
    
    
    n = 3;
     
    model.xs=xs;
    model.ys=ys;
    model.xt=xt;
    model.yt=yt;
    model.xobs=xobs;
    model.yobs=yobs;
    model.robs=robs;
    model.n=n;
    model.xmin=xmin;
    model.xmax=xmax;
    model.ymin=ymin;
    model.ymax=ymax;
    
end