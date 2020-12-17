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

function model=CreateModel(obs, verbose)
    % Domain
    xmin = -40;
    ymin = -40;
    xmax = 40;
    ymax = 40;

    % Start
    xs = -35;
    ys = 35;

    % Target (Destination)
    xt = 35;
    yt = -35;

    % ---------
    % Obstacles
    % ---------

    if obs.load
      xobs = [];
      yobs = [];
      robs = [];

      data = importdata(obs.filename);
      for i = 1:length(data.textdata)
        obs = data.data(i,:);
        assetdata = importdata(sprintf('Data/Obstacles/Assets/%s.txt', data.textdata{i}));
        assetdata = assetdata * obs(5);

        center = [mean(assetdata(:, 1)), mean(assetdata(:, 2))];
        temp = [assetdata ; center];
        dists = squareform(pdist(temp));
        dists = dists(end,:);
        radius = max(dists) + 1;

        obs = data.data(i,:);
        center(1) = center(1) + obs(1);
        center(2) = center(2) + obs(2);

        % theta=linspace(0,2*pi,100);
        % fill(center(1) + radius * cos(theta), center(2) + radius * sin(theta), [0.5 0.7 0.8]);
        % hold on;

        R = [cos(obs(4)) -sin(obs(4)); sin(obs(4)) cos(obs(4))];
        assetdata=assetdata*R;
        % plot(obs(1) + assetdata(:,1), obs(2) + assetdata(:,2))

        xobs = [xobs, center(1)];
        yobs = [yobs, center(2)];
        robs = [robs, radius];
      end
    else

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
    end

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
