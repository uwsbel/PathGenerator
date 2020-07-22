function [] = generateMetrics(S)
% Calculates the metrics for each simulation result

% Input parameter settings parsing
if nargin == 0
    S = struct();
end

if ~isfield(S,'plot') || isempty(S.plot)
    % Don't plot the results
    S.plot = false;
end

if ~isfield(S,'save') || isempty(S.save)
    % Don't save the results
    S.save = false;
end

if ~isfield(S,'verbose') || isempty(S.verbose)
    % Set verbose to off
    S.verbose = false;
end

if ~isfield(S,'results_folder') || isempty(S.results_folder)
    % Set the default folder where results are stored
    S.results_folder = 'Data/Results/';
end

if ~isfield(S,'config') || isempty(S.config)
    % Set the default config to 1F_3L
    S.config = '1L_3F/';
end

if ~isfield(S,'terrain') || isempty(S.terrain)
    % Set the default terrain type to rigid
    S.terrain = 'rigid';
end

if S.save && (~isfield(S,'output_file') || isempty(S.output_file))
    % Set the default output file to be Data/metrics.txt
    S.output_file = 'Data/metrics.txt';
end

if ~isfield(S, 'leaders') || isempty(S.leaders)
   % Leader is only rank 4, by default
   % 4 total agent ranks, where 4 is the leader and 1,2,3 are the followers
   S.leaders = 4;
end

if ~isfield(S, 'followers') || isempty(S.followers)
   % Calculate metrics for all followers by default
   % 4 total agent ranks, where 4 is the leader and 1,2,3 are the followers
   S.followers = [1,2,3];
end

if ~isfield(S, 'max_it') || isempty(S.max_it)
   % Set the maximum number of iterations to 128
   S.max_it = 128;
end

results_dir = sprintf('%s%s%s', S.results_folder, S.config, S.terrain);
if S.verbose
    fprintf('Generating metrics from results stored in %s.\n', results_dir);
end

% output file to save metric data to
ofile = sprintf('Data/metrics_%s.txt', S.terrain);
if S.save
    % delete current file
    fopen(ofile, 'w');
    fprintf(ofile, '#path_idx, rank, path_avg, path_max, speed_avg, speed_max, rel_speed_avg, rel_speed_max, path_length');
end
ofile = fopen(ofile, 'a');

i = 0;
while i < S.max_it
    i = i + 1;

%     if i == 65
%       fprintf('Skipping Bad File...\n');
%       continue;
%     end

    fprintf('Reading %s%d.\n', S.terrain, i);

    while ~isfolder(sprintf('%s/%s%d', results_dir, S.terrain, i))
        i = i + 1;
    end

    % Leader will always be 4 in this case. Should really be a loop
    leader = S.leaders;
    for follower = S.followers
        % Generate leader and follower file paths
        lfile = sprintf('%s/%s%d/output%d.txt', results_dir, S.terrain, i, leader);
        ffile = sprintf('%s/%s%d/output%d.txt', results_dir, S.terrain, i, follower);
        if S.verbose
            fprintf('Leader file: %s.\n', lfile);
            fprintf('Follower file: %s.\n', ffile);
        end

        % calculate metrics
        try
          [mPath, mSpeed, mSpeedRel, path_length] = metrics(lfile, ffile, S.plot);
        catch
          fprintf('Skipping %d.', i);
          continue;
        end

        if S.verbose
          fprintf('File %s%d.\n', S.terrain, i);
          fprintf('Path:\n\tAvg: %f\n\tMax: %f\n\tLength: %f.\n', mPath.avg, mPath.max, path_length);
          fprintf('Speed:\n\tAvg: %f\n\tMax: %f.\n\n', mSpeed.avg, mSpeed.max);
        end

        % save metrics
        if S.save
            fprintf(ofile, '%d,%d,%f,%f,%f,%f,%f,%f,%f\n', i, follower, mPath.avg, mPath.max, mSpeed.avg, mSpeed.max, mSpeedRel.avg, mSpeedRel.max, path_length);
        end
    end
end

fclose(ofile);
