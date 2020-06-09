function [circData, circHandles, frame, S] = bubblebath(S)
% Create a 2D plot of circles (or 'bubbles') with random centers and varying radii.
% Optional Inputs
%   * S is a structure with the following fields, all of which are optional. Their
%     default values and detailed descriptions are within the first block of code.
%     axisHandle: axis handle; an axis is created in a new figure if missing.
%      frameSize: 1x2 vector defining the [x,y] axis lengths (centered at 0).       default [50,50]
%       circSize: 1x2 vector defining [min, max] range of radii.                    default [0.2,5]
%                 1xn vector defining all possible radii when nSizes is NaN.
%         nSizes: Numeric scalar, number of radii to include within circSize.       default 25
%                 NaN: indicates that circSize defines all radii.
% maxCircsPerRad: Maximum number of circles per radius to avoid memory problems.    default 5000
%                 Positive integer or inf.
%     circPoints: Scalar integer greater or equal to 3 describing the number of     default 628
%                 coordinates used to generate each polygon between [0, 2*pi].
%                 Use S.circPoints=3 to draw triangles, = 4 for squares, etc.
%          maxIt: Numeric scalar, max number of attempts to draw circle per radius. default 200
%       edgeType: 0:circle edges can expand outside the frame.                      default 2
%                 1:circles must be entirely inside the frame.
%                 2:circle edges that would expand beyond the frame are cut off.
%        density: Density of circles; numeric scalar [0 < density <= 1].            default 0.7
%        overlap: True: circles can overlap any amount.                             default false
%                 False|0: circle edges can touch but not overlap.
%                 Positive scalar: sets the minimum edge distance (see overlapType).
%                 Negative scalar: sets the maximum edge overlap (see overlapType).
%    overlapType: 'relative': overlap amount is r*overlap, r is radius of cirlce.   default 'relative'
%                 'absolute': overlap amount is constant for all radii.
%      initCircs: 3xn vector defining [x,y,rad] of starter circles                  default NaN
% supressWarning: true: supress internal warnings; false: show warning.             default false
%      drawFrame: true: draw a rectangle showing the frame; false: don't draw it.   default true
%      drawScene: true: draw a the scenes                                           default true
%
% Outputs
%   * circData: [m x 3] matrix of m circles data showing [xCenter, yCenter, radius].
%   * circHandles: [m x 1] vector of handles for each line object / circle.
%   * frame: handle to 'frame' rectangle (GraphicsPlaceHolder if drawFrame is false).
%   * S: A structure with fields listing all parameters used to reproduce the figure
%       also including S.rng which is the random number generator state you can use
%       to reproduce a figure.
%
% Examples
%   Example 1: draw bubblebath plot using default inputs.
%       BUBBLEBATH()
%   Example 2: specify some of the inputs
%       S.circSize = [0.5, logspace(0, 1, 5)];
%       S.nSizes = NaN;
%       S.maxIt = 50;
%       BUBBLEBATH(S)
%   Example 3: Reproduce an exact replica
%       [~,~,~,Sout] = bubblebath();
%       newFig = figure();
%       Sout.axisHandle = axes(newFig);
%       rng(Sout.rng) % seed the random number generator
%       bubblebath(Sout)
%
% Examples of how to use outputs
%   * Calculate the area for each circle: circArea = pi * circData(:,3).^2;
%   * Change color and width of circle edges:  set(circHandles, 'color', 'r', 'LineWidth', 3)
%   * Fill circles with color:
%         c = jet(size(circData,1)); %choose colors
%         for i = 1:size(circData,1)
%             rectangle('Position', [circData(i,1:2)-circData(i,3), circData(i,[3,3])*2], 'Curvature',[1,1], 'FaceColor',c(i,:))
%         end
%   * Count number of circles for each circle size:
%       [g, radius] = findgroups(circData(:,3));
%       count = splitapply(@length, circData(:,3), g);
%       table(radius, count)
%
% Details
%   * For each radius, the number of circles we attempt to draw is a proportion between
%       the circle's area and the frame's area. For 'n' circle sizes, the total area of
%       the frame is evenly split into 'n' portions and the number of circles per circle
%       size is merely framePortion / circleArea, rounded.
%   * Circles are drawn in order of circle size in descending order. The distance
%       between circles is calculated by the distance between all circle centers minus
%       the radii.  That gives the distance between all circle edges.
%   * Since circle centers are chosen randomly, we check that the edges won't overlap
%       with existing circles.  If they do, we choose another set of random centers
%       and that process repeats in a while-loop until either all circle centers
%       do not overlap or we meet a maximum iteration (set by user).  A warning
%       message appears if the while loop expired without drawing all expected
%       circles.
%
% Requires Matlab r2016a or later.
% Copyright (c) 2019, Adam Danz
% All rights reserved
% source: https://www.mathworks.com/matlabcentral/fileexchange/70348-draw-randomly-centered-circles-of-various-sizes
% Contact adam.danz@gmail.com for questions, bugs, suggestions, and high-fives.

% This function was written for a Matlab forum question:
% https://www.mathworks.com/matlabcentral/answers/446114-non-overlapping-random-circles

% Revision history
% vs 1.0.0  02/21/2019  Initial FEX upload
% vs 1.2.0  09/02/2019  Added waitbar, adapted to work with r2016b.
% vs 2.0.0  02/02/2020  Converted to proper function with input options; added several new
%                       options: supressWarning, overlap, overlapType, maxCircsPerRadius,
%                       axisHandle; names of some variables; input validation; documentation.
% vs 2.1.0  02/03/2020  S.nSizes=nan was being overwritten which affected the S output. Fixed.
%                       Updated waitbar message in include max number of circles being drawn.
% vs 2.2.0  04/07/2020  S.circPoints added (ability to draw n-sided polygons). Unrecognized
%                       fields of S will now throw warning rather than an error. fs2 renamed
%                       intlFrame and changed to 2xn mat so S.edgeType=1 fits more bubbles.
%                       Max iterations warning now without backtrace. maxCircsPerRadiuschanged
%                       to maxCircsPerRad but both are accepted for backward compat.

%% Default inputs
% Fill in missing params with default values.
if nargin == 0
    S = struct();
end
if ~isfield(S,'drawScene') || isempty(S.drawScene)
    % Draw the entire scene
    S.drawScene = true;
else
    S.drawFrame = false;
    S.axisHandle = false;
end
if ~isfield(S,'axisHandle') || isempty(S.axisHandle)
    % Create figure & axes
    fh = figure();
    S.axisHandle = axes(fh); % r2016a or later
end
if ~isfield(S,'frameSize') || isempty(S.frameSize)
    % Frame size (arbitrary units): size of axes, centered at (0,0).
    S.frameSize = [50, 50]; %[width, height]
end
if ~isfield(S,'circSize') || isempty(S.circSize)
    % Circle sizes: Select the minimum and maximum radii
    % Radii are linearly spaced. ie: linspace(S.circSize(1), S.circSize(2), S.nSizes).
    % Alternatively, specify all radii directly in a 1xn vector (nSizes must be NaN).
    S.circSize = [.2, 5]; %[smallest, largest] radius
end
if ~isfield(S,'nSizes') || isempty(S.nSizes)
    %number of circle sizes between S.circSize(1) and S.circSize(2)
    % When NaN, circSize specifies radii directly.
    S.nSizes = 25;
end
if isfield(S, 'maxCircsPerRadius')
    % vs 2.0 mistakenly used maxCircsPerRad in help section but the
    % true field name was maxCircsPerRadius.  From vs 2.2 we'll
    % accept both and will change fieldname here.
    S.maxCircsPerRad = S.maxCircsPerRadius;
    S = rmfield(S, 'maxCircsPerRadius');
end
if ~isfield(S,'maxCircsPerRad') || isempty(S.maxCircsPerRad)
    % Set the maximum number of circles per radius.  The squareform
    % function may cause memory problems as this number gets larger
    % (eg ~12000 for 64bit Windows 10).  Use inf for unlimited.
    S.maxCircsPerRad = 5000;
end
if ~isfield(S,'circPoints') || isempty(S.circPoints)
    % Set the number of coordinates used to draw each circle [0,2*pi].
    S.circPoints = 628;
end
if ~isfield(S,'maxIt') || isempty(S.maxIt)
    % max iterations: how many attempts should be made to find circle
    % locations that do not overlap with other circles?
    S.maxIt = 200;
end
if ~isfield(S,'edgeType') || isempty(S.edgeType)
    % Decide whether circles shoule be entirely inside of the frame
    % 0 = Cirlces edges can expand outside of the frame
    % 1 = Cirlces must be entirely inside of the frame
    % 2 = Circle edges that extend beyond the frame are cut off
    S.edgeType = 2;
end
if ~isfield(S,'density') || isempty(S.density)
    % Density of circles:
    % a value greater than 0 but less than or equal to 1.
    % As the value approaches 0, less circles will be drawn.
    S.density = .7;
end
if ~isfield(S,'overlap') || isempty(S.overlap)
    % Allow overlap: True/False or a numeric scalar.
    % When true, bubbles can overlap any ammount.
    % When false or 0, bubble edges can touch but cannot overlap.
    % When numeric & S.overlapType is 'absolute', then all bubble edges will
    % have a minimum edge distance of S.overlap. When S.overlabType
    % is 'relative', then all bubble edges will have a min edge distance
    % of r*S.overlap where r is the radius of the smaller circle
    % being added. Negative values specify the amount of overlap allowed.
    S.overlap = false;
end
if ~isfield(S,'overlapType') || isempty(S.overlapType)
    % See description for overlap.  This parameter is
    % ignored when overlap is 0, true, or false.
    S.overlapType = 'relative';
end
if ~isfield(S,'supressWarning') || isempty(S.supressWarning)
    % When true, the internal warning will be supressed.
    S.supressWarning = false;
end
if ~isfield(S,'drawFrame') || isempty(S.drawFrame)
    % Draw the defined frame around the bubble plot.
    S.drawFrame = true;
end
if ~isfield(S,'initCircs') || isempty(S.initCircs)
    % No predefined circles
    S.initCircs = NaN;
end

% input validation
if S.axisHandle ~= false
    assert(ishghandle(S.axisHandle,'axes'),'axisHandle must be an axis handle.')
end
validateattributes(S.frameSize,{'double'},{'size',[1,2],'>',0},mfilename,'frameSize')
assert(isnumeric(S.nSizes) && isscalar(S.nSizes) && (isnan(S.nSizes) || mod(S.nSizes,1)==0),...
    'nSizes is expected to be a scalar, nonzero, positive integer or NaN.')
if isnan(S.nSizes)
    assert(isrow(S.circSize) && isnumeric(S.circSize) && all(~isinf(S.circSize)) && all(S.circSize>0), ...
        'When nSizes is NaN, circSize is expected to be a 1xn numeric vector, with finite, positive, nonzero values.')
elseif isnumeric(S.circSize) && numel(S.circSize) > 2
    error('When defining all possible radii in S.circSize, S.nSizes must be NaN.')
else
    validateattributes(S.circSize,{'double'},{'size',[1,2],'>',0,'increasing'},mfilename,'circSize')
end
assert(isnumeric(S.maxCircsPerRad) && isscalar(S.maxCircsPerRad) && S.maxCircsPerRad>0 ...
    && (mod(S.maxCircsPerRad,1)==0 || isinf(S.maxCircsPerRad)), 'maxCircsPerRad must be a positive, non-zero integer or inf.')
validateattributes(S.circPoints,{'double'},{'numel',1,'>=',3,'integer'},mfilename,'circPoints')
validateattributes(S.maxIt,{'double'},{'numel',1,'>',0,'integer'},mfilename,'maxIt')
assert(ismember(S.edgeType,0:2),'edgeType must be 0, 1, or 2.')
validateattributes(S.density,{'double'},{'numel',1,'>',0,'<=',1},mfilename,'density')
assert((islogical(S.overlap)||isnumeric(S.overlap))&&numel(S.overlap),'overlap must be true, false, or a numeric scalar value.')
assert(any(strcmpi(S.overlapType,{'absolute','relative'})),'overlapType must be either ''relative'' or ''absolute''.')
validateattributes(S.supressWarning,{'logical'},{'numel',1},mfilename,'supressWarning')
validateattributes(S.drawFrame,{'logical'},{'numel',1},mfilename,'drawFrame')

% Check for included parameters (fields of S) that are not recognized.
acceptedFields = {'axisHandle','frameSize','circSize','nSizes','circPoints','maxIt','edgeType','density',...
    'overlap','overlapType','supressWarning','drawFrame','maxCircsPerRad'};
params = fields(S);
foreignFields = ismember(params,acceptedFields);
if ~all(foreignFields) && ~S.supressWarning
	warning('Parameter(s) [%s] not recognized and will be ignored.', strjoin(params(~foreignFields),', '))
end

% Store rng state for reproducibility
S.rng = rng();

%% Add soap
if S.drawScene
    hold(S.axisHandle,'on')
    axis(S.axisHandle,'equal') %set aspect ratio to 1:1
    xlim(S.axisHandle, S.frameSize(1)/2 * [-1.05, 1.05]) %with some extra space
    ylim(S.axisHandle, S.frameSize(2)/2 * [-1.05, 1.05]) %with some extra space
end

% determine minimum distance between circles allowed
if islogical(S.overlap) && ~S.overlap
    minDist = 0;
elseif islogical(S.overlap) && S.overlap
    minDist = -inf;
else
    minDist = S.overlap;
end

% determine circle sizes (largest to smallest)
if isnan(S.nSizes)
    nSizes = numel(S.circSize);
    r = sort(S.circSize,'descend');
else
    nSizes = S.nSizes;
    r = linspace(S.circSize(2), S.circSize(1), nSizes);
end

% Identify the internal frame of possible circle centers (intlFrame)
switch S.edgeType
    case 0  % Cirlce edges can expand outside of the frame
        intlFrame = repmat(S.frameSize(:), 1, nSizes);
        clearBorder = false;
    case 1  % Circle must be entirely inside of the frame
        intlFrame = bsxfun(@minus, S.frameSize(:), r*2);
        clearBorder = false;
    case 2  % Circle edges end at the frame
        intlFrame = repmat(S.frameSize(:), 1, nSizes);
        clearBorder = true;
end

% adjust min distance between circles based on overlap type
if strcmpi(S.overlapType,'absolute')
    minDist = repmat(minDist,size(r));
elseif strcmpi(S.overlapType,'relative')
    minDist = r*minDist;
end

% determine approximate number of circles to draw for each size
frameArea = prod(S.frameSize);
circAreas = pi * r.^2;
d = ceil(ceil((frameArea / nSizes) ./ circAreas) * S.density); %Old: (2:nSizes+1).^2;
if any(d > S.maxCircsPerRad) && ~S.supressWarning
    warning('The maximum number of circles for radii [%s] have been limited to %d.', ...
        regexprep(strtrim(num2str(r(d>S.maxCircsPerRad))),' +',' '), S.maxCircsPerRad)
end
d(d > S.maxCircsPerRad) = S.maxCircsPerRad;

% Throw error if largest circle is larger than frame
assert(max(circAreas) <= frameArea, 'The area of the largest circle (%.1f) is larger than the frame''s area (%.1f)', max(circAreas), frameArea)

% Loop through each circle size
if ~isnan(S.initCircs)
  circdata = S.initCircs;
else
  circdata = []; %[xCenter, yCenter, radius] of each drawn circle
end
h = cell(nSizes,1); % handles to the line objects for each circle
wb = waitbar(0,'initializing...','name',mfilename);
originalWarnState = warning('backtrace');
warning backtrace off
for i = 1:nSizes
    % Reset circle & iteration count
    cCount = 0;
    iCount = 0;
    xRand = zeros(d(i),1);
    yRand = xRand;
    isOK = false(size(xRand));
    % Keep drawing random coordinates until either all circs
    % are drawn or we reach max number of attempts
    while cCount < d(i) && iCount < S.maxIt
        % Randomly choose center coordinates
        xRand(~isOK) = (rand(d(i)-cCount,1) - 0.5) * intlFrame(1,i);
        yRand(~isOK) = (rand(d(i)-cCount,1) - 0.5) * intlFrame(2,i);
        % determine if new circles overlap with others or themselves
        xyr = [circdata; [xRand, yRand, repmat(r(i), size(yRand))]];
        radMat = repmat(xyr(:,3),1,size(xyr,1)) + repmat(xyr(:,3).',size(xyr,1),1); %changed 190901 to work with r2016a
        dist = tril(squareform(pdist(xyr(:,1:2))) - radMat, -1);
        if isempty(dist)
            isOK = true; %when xyr only has 1 row
        else
            dist(triu(true(size(dist)))) = inf;
            isOK = all(dist(size(circdata,1)+1:end,:) >= minDist(i), 2);
        end
        cCount = sum(isOK);  %cirlce count for current radius
        iCount = iCount + 1; %iteration count
        % Update waitbar
        if ishghandle(wb)
            waitbar(max(iCount/S.maxIt,cCount >= d(i)),wb,sprintf('Trying to find space for up to %d circles with radius = %.2f\nFinal radius: %.2f',d(i),r(i),r(end)));
        end
    end
    % If we had to quit searching, throw warning.
    if iCount >= S.maxIt && ~S.supressWarning
        warning('Max iteration reached. %d of the requested %d circles drawn for radius %.3f', cCount, d(i), r(i))
    end
    % Store all final circle data
    circdata = [circdata; [xRand(isOK), yRand(isOK), repmat(r(i), sum(isOK), 1)]]; %#ok<AGROW>
    % Draw circles
    if any(isOK)
        if S.drawScene
            h{i} = drawcircles([xRand(isOK), yRand(isOK), repmat(r(i), sum(isOK))], clearBorder, S);
        end
    end

end
warning(originalWarnState) %return original warning state
% Draw frame
if S.drawFrame
    frame = rectangle(S.axisHandle, 'position', [-S.frameSize/2, S.frameSize], 'LineWidth', 2);
else
    frame = gobjects(1);
end
circHandles = [h{:}]';

% Remove waitbar
if ishghandle(wb)
    delete(wb)
end

if nargout > 0
    % Produce output only when requested
    circData = circdata;
end

function h = drawcircles(xyr, clearBorder, S)
% Draw circle given center and radius
ang = linspace(0, 2*pi, S.circPoints+1);
xp = xyr(:,3)*cos(ang) + repmat(xyr(:,1),1,numel(ang)); %changed 190901 to work with r2016a
yp = xyr(:,3)*sin(ang) + repmat(xyr(:,2),1,numel(ang)); %changed 190901 to work with r2016a
if clearBorder
    % remove data outside of frame, if requested
    xp(abs(xp) > S.frameSize(1)/2) = NaN;
    yp(abs(yp) > S.frameSize(2)/2) = NaN;
end

h = plot(S.axisHandle, xp', yp', 'k')';
