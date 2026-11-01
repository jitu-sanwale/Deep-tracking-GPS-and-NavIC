function trajTable = generateFig7Trajectory(outFile, maneuverStart, jerkPeak_gps, V0, dt)
%GENERATEFIG7TRAJECTORY Build a 60 s single continuous ground-truth
%   trajectory CSV containing one high-dynamics maneuver that reproduces
%   the acceleration/jerk *shape* of Won, Pany & Eissfeller, "Iterative
%   Maximum Likelihood Estimators for High-Dynamic GNSS Signal Tracking,"
%   IEEE Trans. Aerospace and Electronic Systems 48(4), Oct 2012, Fig. 7
%   ("Signal dynamic truth under condition of high dynamics"), scaled to
%   a harder peak jerk than the paper's 100 g/s (default here: 120 g/s).
%
%   Fig. 7's acceleration profile (read directly off the published
%   figure, not just its caption text) is nine piecewise-constant-jerk
%   segments spanning a fixed 7.0 s:
%     [0.0,2.0]s  cruise            (init lead-in, a = 0)
%     [2.0,2.5]s  ramp,  jerk=+Jp   (a: 0 -> +Ap)
%     [2.5,3.0]s  hold              (a = +Ap)
%     [3.0,4.0]s  ramp,  jerk=-Jp   (a: +Ap -> -Ap, crosses 0 at 3.5s)
%     [4.0,4.5]s  hold              (a = -Ap)
%     [4.5,5.5]s  ramp,  jerk=+Jp   (a: -Ap -> +Ap, crosses 0 at 5.0s)
%     [5.5,6.0]s  hold              (a = +Ap)
%     [6.0,6.5]s  ramp,  jerk=-Jp   (a: +Ap -> 0)
%     [6.5,7.0]s  cruise            (tail, a = 0)
%   where Jp = jerkPeak_gps * 9.81 [m/s^3] and Ap = Jp * 0.5s (the fixed
%   0.5 s ramp/hold timing from the paper is preserved unchanged; only
%   the jerk/accel magnitude is scaled). At the paper's own 100 g/s this
%   reproduces Ap = 50g exactly, matching Fig. 7; at the default 120 g/s
%   used here it reaches Ap = 60g.
%
%   This maneuver is embedded, unmodified in shape, at
%   [maneuverStart, maneuverStart+7] inside an otherwise flat 60 s
%   trajectory: constant-velocity cruise from 0 to maneuverStart, the 7 s
%   maneuver, then constant-velocity cruise (at whatever velocity the
%   maneuver left it at - the +Ap/-Ap segments are not perfectly
%   symmetric in duration, exactly as in the real Fig. 7) out to 60 s.
%
%   Unlike generateHighDynamicsTrajectory.m (independent jerk-doublets
%   that always return acceleration to exactly 0 at every segment
%   boundary) and generateLoseLockTrajectory.m (a single ramp-hold-ramp
%   pulse), this generator carries acceleration state - not just
%   position/velocity - across segment boundaries, since three of Fig.
%   7's hold segments sit at a nonzero acceleration and two of its ramps
%   cross zero acceleration mid-segment. Position/velocity/acceleration
%   are still evaluated from EXACT closed-form per-segment polynomials
%   (not numerical integration of the discontinuous jerk), so there is no
%   accumulated drift regardless of the sample interval dt.
%
%   trajTable = generateFig7Trajectory(outFile, maneuverStart, jerkPeak_gps, V0, dt)
%
%   Inputs (all optional):
%       outFile        - CSV path to write (default 'Fig7Trajectory.csv')
%       maneuverStart  - [s] when the 7 s maneuver begins (default 20,
%                        i.e. the maneuver spans [20, 27] s)
%       jerkPeak_gps   - peak jerk magnitude [g/s] (default 120; the
%                        paper's own scenario is 100)
%       V0             - baseline East ground speed [m/s] before the
%                        maneuver (default 50, matches
%                        generateHighDynamicsTrajectory.m's convention)
%       dt             - trajectory sample interval [s] (default 0.02 s,
%                        i.e. 50 Hz)
%
%   Output:
%       trajTable      - the same [time,lat,lon,height,roll,pitch,heading]
%                        matrix written to outFile, for inspection/plotting.

if nargin < 1 || isempty(outFile),       outFile       = 'Fig7Trajectory.csv'; end
if nargin < 2 || isempty(maneuverStart), maneuverStart = 20;                   end
if nargin < 3 || isempty(jerkPeak_gps),  jerkPeak_gps  = 120;                  end
if nargin < 4 || isempty(V0),            V0            = 50;                  end
if nargin < 5 || isempty(dt),            dt            = 0.02;                end

TOTAL_DUR    = 60;  % s, fixed total trajectory length
MANEUVER_DUR = 7;   % s, fixed - matches the paper's own Fig. 7 time axis
if maneuverStart < 0 || maneuverStart + MANEUVER_DUR > TOTAL_DUR
    error(['generateFig7Trajectory: maneuverStart=%.2f places the 7 s ' ...
        'maneuver outside the [0, %d] s trajectory.'], maneuverStart, TOTAL_DUR);
end

gAccel = 9.81;
Jp = jerkPeak_gps * gAccel;   % m/s^3, peak jerk magnitude
Ap = Jp * 0.5;                % m/s^2, peak accel reached after a 0.5 s ramp

%% Build the segment list: {duration, jerk} ---------------------------------
tailCruiseDur = TOTAL_DUR - maneuverStart - MANEUVER_DUR;
segs = struct('duration', {}, 'jerk', {});
segs(end+1) = struct('duration', maneuverStart, 'jerk', 0);  % lead-in cruise
segs(end+1) = struct('duration', 2.0,           'jerk', 0);  % Fig.7's own init
segs(end+1) = struct('duration', 0.5,           'jerk', +Jp); % ramp up
segs(end+1) = struct('duration', 0.5,           'jerk', 0);   % hold at +Ap
segs(end+1) = struct('duration', 1.0,           'jerk', -Jp); % ramp thru 0 to -Ap
segs(end+1) = struct('duration', 0.5,           'jerk', 0);   % hold at -Ap
segs(end+1) = struct('duration', 1.0,           'jerk', +Jp); % ramp thru 0 to +Ap
segs(end+1) = struct('duration', 0.5,           'jerk', 0);   % hold at +Ap
segs(end+1) = struct('duration', 0.5,           'jerk', -Jp); % ramp down to 0
segs(end+1) = struct('duration', 0.5,           'jerk', 0);   % Fig.7's own tail (completes its 7.0s window)
segs(end+1) = struct('duration', tailCruiseDur, 'jerk', 0);   % outer cruise back to 60s total

segStart = zeros(1, numel(segs));
for i = 2:numel(segs)
    segStart(i) = segStart(i-1) + segs(i-1).duration;
end
totalDur = segStart(end) + segs(end).duration;  % should equal TOTAL_DUR

%% Walk the segments, carrying exact state (x,v,a) across boundaries --------
xAtStart = zeros(1, numel(segs));
vAtStart = zeros(1, numel(segs));
aAtStart = zeros(1, numel(segs));
vAtStart(1) = V0;
maneuverLog = {};
for i = 1:numel(segs)
    x0 = xAtStart(i); v0 = vAtStart(i); a0 = aAtStart(i);
    J = segs(i).jerk; D = segs(i).duration;
    aEnd = a0 + J*D;
    vEnd = v0 + a0*D + 0.5*J*D^2;
    xEnd = x0 + v0*D + 0.5*a0*D^2 + (J/6)*D^3;
    if J ~= 0 || a0 ~= 0
        maneuverLog{end+1} = sprintf(...
            ['  t=%6.2f-%6.2fs: jerk=%+9.1f m/s^3, accel %+7.1f -> %+7.1f m/s^2 ' ...
             '(%+5.1fg -> %+5.1fg), speed %8.2f -> %8.2f m/s'], ...
            segStart(i), segStart(i)+D, J, a0, aEnd, a0/gAccel, aEnd/gAccel, v0, vEnd); %#ok<AGROW>
    end
    if i < numel(segs)
        xAtStart(i+1) = xEnd; vAtStart(i+1) = vEnd; aAtStart(i+1) = aEnd;
    end
end

%% Sample each segment on the fine time grid ---------------------------------
t = 0:dt:totalDur;
if t(end) < totalDur - eps(totalDur)*10
    t = [t, totalDur]; %#ok<AGROW>
end
posE = zeros(size(t));
vel  = zeros(size(t));
acc  = zeros(size(t)); %#ok<NASGU> % kept for callers that want the accel history too

for i = 1:numel(segs)
    tSegStart = segStart(i);
    tSegEnd   = tSegStart + segs(i).duration;
    if i < numel(segs)
        inSeg = t >= tSegStart & t < tSegEnd;
    else
        inSeg = t >= tSegStart & t <= tSegEnd + eps;
    end
    tRel = t(inSeg) - tSegStart;
    x0 = xAtStart(i); v0 = vAtStart(i); a0 = aAtStart(i); J = segs(i).jerk;
    posE(inSeg) = x0 + v0*tRel + 0.5*a0*tRel.^2 + (J/6)*tRel.^3;
    vel(inSeg)  = v0 + a0*tRel + 0.5*J*tRel.^2;
    acc(inSeg)  = a0 + J*tRel;
end

%% Convert East offset to geodetic lat/lon (flat-Earth tangent plane) ------
% Same origin/altitude convention as generateHighDynamicsTrajectory.m /
% generateLoseLockTrajectory.m, so this trajectory drops into the same
% satellite-visibility geometry as the project's other stress tests.
lat0 = 20.1085; lon0 = 73.9532; h0 = 1000;
Rearth = 6371000; % m, mean Earth radius
lat = lat0 * ones(size(t));
lon = lon0 + (posE / (Rearth * cosd(lat0))) * (180/pi);
height = h0 * ones(size(t));
roll = zeros(size(t));
pitch = zeros(size(t));
heading = 90 * ones(size(t));  % due East throughout

trajTable = [t(:), lat(:), lon(:), height(:), roll(:), pitch(:), heading(:)];

%% Write CSV ------------------------------------------------------------------
fid = fopen(outFile, 'w');
if fid < 0
    error('generateFig7Trajectory: could not open %s for writing.', outFile);
end
% Header line ends in \r\n to match the line ending writematrix's
% 'WriteMode','append' produces for the data rows below on Windows/MATLAB
% (see the same fix in generateHighDynamicsTrajectory.m).
fprintf(fid, ['Time [sec],Latitude [deg],Longitude [deg],Height [meter],' ...
    'Roll Angle [deg],Pitch Angle [deg],Heading [deg]\r\n']);
fclose(fid);
writematrix(trajTable, outFile, 'WriteMode', 'append');

fprintf('Fig.7-style high-dynamics trajectory written: %s (%d points, %.1f s total)\n', ...
    outFile, size(trajTable,1), totalDur);
fprintf('Maneuver window: %.1f-%.1fs, peak jerk +/-%.1f m/s^3 (%.0f g/s), peak accel +/-%.1f m/s^2 (%.0fg)\n', ...
    maneuverStart, maneuverStart+MANEUVER_DUR, Jp, jerkPeak_gps, Ap, Ap/gAccel);
for i = 1:numel(maneuverLog)
    fprintf('%s\n', maneuverLog{i});
end
fprintf('  Speed: %.2f m/s (t=0) -> %.2f m/s (post-maneuver cruise, t=%.1f-%.1fs)\n', ...
    vel(1), vel(end), maneuverStart+MANEUVER_DUR, totalDur);
