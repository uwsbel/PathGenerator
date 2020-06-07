
num_trials = 5;
verbose = false;
render = true;

num_bezier_points = 30;

% ----------------

addpath(genpath('PathPlanning/'));
addpath(genpath('SampleCurve/'));
addpath(genpath('Bezier/'));

% ----------------

boolString = {'false', 'true'};
idx = 0;
for it = 1:num_trials
    
    % Generate optimal path for random obstacles
    fprintf('Trial %d... ', it);
    tstart = tic;
    [sol, model] = pso(verbose);
    feasible = sol.IsFeasible;
    points = [sol.XS' sol.YS'];
    obst = [model.xobs' model.yobs' model.robs'];
    fprintf('feasible? %s   [time: %.2f]\n', boolString{feasible+1}, toc(tstart));
    
    % TODO:  a more precise feasability test?
    
    if ~feasible
        continue;
    end

    % Create Bezier curve 
    p = interparc(num_bezier_points, sol.xx, sol.yy);
    data = generateBezier([p zeros(size(p,1),1)]);
    
    % Render obstacles, path, and Bezier curve
    if render
        figure;
        hold on
        
        PlotSolution(sol, model);
        renderBezier(data, 100);
        
        axis equal;
        grid on
        box on
        xlim([-50 50])
        ylim([-50 50])
        
        drawnow
    end
    
    % Write files for this trial
    idx = idx + 1;
    
    % Path points
    fp = fopen(sprintf('points_%d.dat', idx), 'w');
    for i = 1:size(points,1)
        fprintf(fp, '%f %f\n', points(i,1), points(i,2));
    end
    fclose(fp);
    
    % Obstacle information
    fo = fopen(sprintf('obst_%d.dat', idx), 'w');
    for i = 1:size(obst,1)
        fprintf(fo, '%f %f %f\n', obst(i,1), obst(i,2), obst(i,3));
    end
    fclose(fo);
    
    % Bezier curve
    writeBezier(data, sprintf('path_%d.dat', idx));

end

fclose('all');

