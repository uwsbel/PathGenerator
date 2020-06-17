
num_trials = 750;
verbose = false;
render = false;
save = true;

num_bezier_points = 30;

% ----------------

addpath(genpath('PathPlanning/'));
addpath(genpath('SampleCurve/'));
addpath(genpath('Bezier/'));
addpath(genpath('Obstacles/'));

% ----------------

boolString = {'false', 'true'};
idx = 1;
for it = 1:num_trials

    % Generate optimal path for random obstacles
    fprintf('Trial %d... ', it);
    tstart = tic;
    [sol, model] = pso(verbose);
    feasible = sol.IsFeasible;
    points = [sol.XS' sol.YS'];
    obst = [model.xobs' model.yobs' model.robs'];
    fprintf('feasible? %s   [time: %.2f]\n', boolString{feasible+1}, toc(tstart));

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

    % TODO:  a more precise feasability test?

    if ~feasible
        continue;
    end

    if save
        % Write files for this trial
        folder = "Data/Paths/";
        while isfile(sprintf('%spath_%d.dat', folder, idx))
            idx = idx + 1;
        end

        % Path points
        fp = fopen(sprintf('%spoints_%d.dat', folder, idx), 'w');
        for i = 1:size(points,1)
            fprintf(fp, '%f %f\n', points(i,1), points(i,2));
        end
        fclose(fp);

        % Obstacle information
        fo = fopen(sprintf('%sobst_%d.dat', folder, idx), 'w');
        for i = 1:size(obst,1)
            fprintf(fo, '%f %f %f\n', obst(i,1), obst(i,2), obst(i,3));
        end
        fclose(fo);

        % Bezier curve
        writeBezier(data, sprintf('%spath_%d.dat', folder, idx));

        fprintf('Figure #%d saved.\n', idx);
    end

end

fclose('all');
