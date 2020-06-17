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

function model=CreateModel(verbose)
    % Domain
    xmin = -46;
    ymin = -46;
    xmax = 46;
    ymax = 46;

    % Start
    xs = -46;
    ys = -46;

    % Target (Destination)
    xt = 46;
    yt = 46;

    % ---------
    % Obstacles
    % ---------

    % Fixed (equally spaced between start and target)
    nfixed = 4;                         %  CAN BE CHANGED
    if nfixed
      foo = linspace(xs, xt, nfixed + 2);
      xfixed = foo(2:end-1);
      foo = linspace(ys, yt, nfixed + 2);
      yfixed = foo(2:end-1);
      rfixed = 8 * ones(1, nfixed);       %  CAN BE CHANGED
    end

    % Initialize obstacle generation settings
    S.circSize = [5,6,7];
    S.frameSize = [85,85];
	S.maxIt = 50;
    S.nSizes = NaN;
    S.supressWarning = ~verbose;
    S.maxCircsPerRad = 5;
    S.drawScene = false;
    S.overlap = -1;
    if nfixed
      S.initCircs = [xfixed; yfixed; rfixed]';
    end

    % Create the obstacles
    circData = bubblebath(S);

    % All obstacles
    xobs = circData(:,1)';
    yobs = circData(:,2)';
    robs = circData(:,3)';

    if verbose
      fprintf('Number of Obstacles: %d.\n', length(xobs))
    end

    % ===================================================================


    n = 5;

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
