function trajTable = generateHighDynamicsTrajectory(outFile, jerkValues, dt)
%GENERATEHIGHDYNAMICSTRAJECTORY Build and write a single continuous
%   ground-truth trajectory CSV consisting of a sequence of jerk-doublet
%   maneuvers of increasing severity, for single-satellite GNSS
%   tracking-loop dynamic-stress testing (jerk = 3rd time-derivative of
%   position, i.e. rate of change of acceleration).
%
%   A "doublet" is the standard flight-test excitation: jerk is held at
%   +J for a half-duration Th, then -J for the next Th, returning
%   acceleration smoothly to zero while leaving a net, permanent change
%   in velocity of J*Th^2 (avoids an unrealistic jerk-discontinuous
%   acceleration step). Motion here is 1D (due East) at constant
%   altitude, so the injected jerk/acceleration/velocity maps directly
%   onto along-track range-rate dynamics as seen by a satellite roughly
%   along that heading - exactly what stresses a PLL/FLL/DLL.
%
%   Position/velocity/acceleration are evaluated from EXACT closed-form
%   per-segment polynomials (not numerical integration of the
%   discontinuous jerk signal), so there is no accumulated drift between
%   segments regardless of the sample interval dt - dt only controls the
%   CSV's time resolution.
%
%   Six segments, one per jerkValues entry (default
%   [0 20 40 80 100 120] m/s^3), each starting at a distinct timestamp
%   along ONE continuous trajectory (matching one IF data file): the
%   first (jerk=0) is a plain constant-velocity baseline with no
%   maneuver; each subsequent entry is an 8 s constant-velocity settling
%   run (so the tracking loop can re-lock/stabilize before the next
%   stress event) followed by a 2 s (Th=1 s) jerk doublet of that
%   magnitude, alternating sign so ground speed stays bounded (net
%   delta-v +20, -40, +80, -100, +120 m/s) rather than growing without
%   bound. An 8 s constant-velocity tail follows the last maneuver.
%
%   trajTable = generateHighDynamicsTrajectory(outFile, jerkValues, dt)
%
%   Inputs (all optional):
%       outFile     - CSV path to write (default 'HighDynamicsTrajectory.csv')
%       jerkValues  - 1xN jerk magnitudes [m/s^3], first is the baseline
%                    (should be 0) (default [0 20 40 80 100 120])
%       dt          - trajectory sample interval [s] (default 0.02 s,
%                    i.e. 50 Hz)
%
%   Output:
%       trajTable   - the same [time,lat,lon,height,roll,pitch,heading]
%                    matrix written to outFile, for inspection/plotting.

if nargin < 1 || isempty(outFile),    outFile    = 'HighDynamicsTrajectory.csv'; end
if nargin < 2 || isempty(jerkValues), jerkValues = [0 20 40 80 100 120];         end
if nargin < 3 || isempty(dt),         dt         = 0.02;                        end

%% Origin and cruise parameters -------------------------------------------
lat0 = 20.1085;         % deg - same origin as static_Trajectory.csv/projectile_Vivek.csv
lon0 = 73.9532;         % deg
h0   = 1000;            % m, constant altitude throughout
V0   = 50;              % m/s, initial (baseline) East ground speed
Th   = 1;               % s, half-duration of each jerk doublet (2 s total)
cruiseDur = 8;          % s, constant-velocity settle time around each maneuver
initialCruiseDur = 10;  % s, baseline segment before the first maneuver

nSeg = numel(jerkValues);
% +1,-1,+1,-1,... for the maneuvers after the baseline (jerkValues(2:end))
signPattern = (-1) .^ (0:nSeg-2);

%% Build the segment list: {type, duration, jerk} --------------------------
segs = struct('type', {}, 'duration', {}, 'jerk', {});
segs(end+1) = struct('type', 'cruise', 'duration', initialCruiseDur, 'jerk', 0);
for k = 2:nSeg
    segs(end+1) = struct('type', 'cruise', 'duration', cruiseDur, 'jerk', 0); %#ok<AGROW>
    segs(end+1) = struct('type', 'doublet', 'duration', 2*Th, ...
        'jerk', signPattern(k-1) * jerkValues(k)); %#ok<AGROW>
end
segs(end+1) = struct('type', 'cruise', 'duration', cruiseDur, 'jerk', 0);

segStart = zeros(1, numel(segs));
for i = 2:numel(segs)
    segStart(i) = segStart(i-1) + segs(i-1).duration;
end
totalDur = segStart(end) + segs(end).duration;

%% Walk the segments, computing exact state (x,v) at each boundary --------
% xAtStart(i), vAtStart(i): East position/velocity at the start of segs(i).
xAtStart = zeros(1, numel(segs));
vAtStart = zeros(1, numel(segs));
vAtStart(1) = V0;
maneuverLog = {};
for i = 1:numel(segs)
    x0 = xAtStart(i);
    v0 = vAtStart(i);
    if strcmp(segs(i).type, 'cruise')
        xEnd = x0 + v0 * segs(i).duration;
        vEnd = v0;
    else % doublet
        J = segs(i).jerk;
        % Phase 1 end (tau=Th): a1=J*Th, v1, x1
        v1 = v0 + 0.5*J*Th^2;
        x1 = x0 + v0*Th + (J/6)*Th^3;
        a1 = J*Th;
        % Phase 2 end (sigma=Th)
        vEnd = v1 + a1*Th - 0.5*J*Th^2;
        xEnd = x1 + v1*Th + 0.5*a1*Th^2 - (J/6)*Th^3;
        maneuverLog{end+1} = sprintf(...
            '  t=%.1f-%.1fs: jerk doublet J=%+g m/s^3 (peak accel %+.1f m/s^2, net delta-v %+.1f m/s, %.1f->%.1f m/s)', ...
            segStart(i), segStart(i)+segs(i).duration, J, a1, J*Th^2, v0, vEnd); %#ok<AGROW>
    end
    if i < numel(segs)
        xAtStart(i+1) = xEnd;
        vAtStart(i+1) = vEnd;
    end
end

%% Sample each segment on the fine time grid -------------------------------
t = 0:dt:totalDur;
if t(end) < totalDur
    t = [t, totalDur]; %#ok<AGROW>
end
posE = zeros(size(t));
vel  = zeros(size(t));

for i = 1:numel(segs)
    tSegStart = segStart(i);
    tSegEnd   = tSegStart + segs(i).duration;
    if i < numel(segs)
        inSeg = t >= tSegStart & t < tSegEnd;
    else
        inSeg = t >= tSegStart & t <= tSegEnd + eps; % include the very last sample
    end
    tRel = t(inSeg) - tSegStart;
    x0 = xAtStart(i);
    v0 = vAtStart(i);
    if strcmp(segs(i).type, 'cruise')
        posE(inSeg) = x0 + v0*tRel;
        vel(inSeg)  = v0;
    else
        J = segs(i).jerk;
        half1 = tRel <= Th;
        % Phase 1 (0<=tau<=Th)
        tau = tRel(half1);
        posE1 = x0 + v0*tau + (J/6)*tau.^3;
        vel1  = v0 + 0.5*J*tau.^2;
        % Phase 2 (Th<tau<=2Th), sigma = tau-Th
        v1 = v0 + 0.5*J*Th^2;
        x1 = x0 + v0*Th + (J/6)*Th^3;
        a1 = J*Th;
        sigma = tRel(~half1) - Th;
        posE2 = x1 + v1*sigma + 0.5*a1*sigma.^2 - (J/6)*sigma.^3;
        vel2  = v1 + a1*sigma - 0.5*J*sigma.^2;

        segPosE = zeros(size(tRel));
        segVel  = zeros(size(tRel));
        segPosE(half1)  = posE1;  segPosE(~half1) = posE2;
        segVel(half1)   = vel1;   segVel(~half1)  = vel2;
        posE(inSeg) = segPosE;
        vel(inSeg)  = segVel;
    end
end

%% Convert East offset to geodetic lat/lon (flat-Earth tangent plane) ------
% Adequate for offsets of this scale (a few km); latitude is unchanged
% since motion is purely East.
Rearth = 6371000; % m, mean Earth radius
lat = lat0 * ones(size(t));
lon = lon0 + (posE / (Rearth * cosd(lat0))) * (180/pi);
height = h0 * ones(size(t));
roll = zeros(size(t));
pitch = zeros(size(t));
heading = 90 * ones(size(t));  % due East throughout (ground speed never reverses sign)

trajTable = [t(:), lat(:), lon(:), height(:), roll(:), pitch(:), heading(:)];

%% Write CSV ----------------------------------------------------------------
fid = fopen(outFile, 'w');
if fid < 0
    error('generateHighDynamicsTrajectory: could not open %s for writing.', outFile);
end
% Header line ends in \r\n (not just \n) to match the line ending that
% writematrix's 'WriteMode','append' produces for the data rows below on
% Windows/MATLAB - otherwise the two writes' line endings mismatch at the
% seam and leave a spurious blank line as the file's first "data" row.
fprintf(fid, ['Time [sec],Latitude [deg],Longitude [deg],Height [meter],' ...
    'Roll Angle [deg],Pitch Angle [deg],Heading [deg]\r\n']);
fclose(fid);
writematrix(trajTable, outFile, 'WriteMode', 'append');

fprintf('High-dynamics trajectory written: %s (%d points, %.1f s, jerk values %s m/s^3)\n', ...
    outFile, size(trajTable,1), totalDur, mat2str(jerkValues));
for i = 1:numel(maneuverLog)
    fprintf('%s\n', maneuverLog{i});
end
fprintf('  Speed: %.1f m/s (start) -> %.1f m/s (end), peak %.1f m/s\n', vel(1), vel(end), max(vel));
