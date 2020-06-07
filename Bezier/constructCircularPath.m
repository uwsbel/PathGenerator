function data = constructCircularPath(run, radius, offset, nturns)
% Construct a Bezier curve to approximate a circular path of given radius
% (centered at x = radius + offset, y = 0, z = 0.1) preceded by a
% straigh-line run of given length.  The path wraps 'nturns' times
% counter-clockwise.
%
% Radu Serban, 9/28/2107

% Approximate circular path using 4 points
factor = radius * (4/3) * tan(pi/8);
z = 0.1;

P1 = [radius + offset, -radius, z];
P1_in = P1 - [factor, 0, 0];
P1_out = P1 + [factor, 0, 0];

P2 = [2 * radius + offset, 0, z];
P2_in = P2 - [0, factor, 0];
P2_out = P2 + [0, factor, 0];

P3 = [radius + offset, radius, z];
P3_in = P3 + [factor, 0, 0];
P3_out = P3 - [factor, 0, 0];

P4 = [offset, 0, z];
P4_in = P4 + [0, factor, 0];
P4_out = P4 - [0, factor, 0];

% Start point
P0 = [radius + offset - run, -radius, z];
P0_in = P0;
P0_out = P0 + [factor, 0, 0];

% Set up output structure
n = 1 + 4 * nturns;
data.n = n;
data.p = zeros(n,3);
data.in = zeros(n,3);
data.out = zeros(n,3);

% Load initial point
data.p(1,:) = P0;
data.in(1,:) = P0_in;
data.out(1,:) = P0_out;

% Load points around circle
for i = 0:nturns-1
   data.p(1 + i*4 + 1, :) = P1;
   data.in(1 + i*4 + 1, :) = P1_in;
   data.out(1 + i*4 + 1, :) = P1_out;
   
   data.p(1 + i*4 + 2, :) = P2;
   data.in(1 + i*4 + 2, :) = P2_in;
   data.out(1 + i*4 + 2, :) = P2_out;
   
   data.p(1 + i*4 + 3, :) = P3;
   data.in(1 + i*4 + 3, :) = P3_in;
   data.out(1 + i*4 + 3, :) = P3_out;
  
   data.p(1 + i*4 + 4, :) = P4;
   data.in(1 + i*4 + 4, :) = P4_in;
   data.out(1 + i*4 + 4, :) = P4_out;
end
