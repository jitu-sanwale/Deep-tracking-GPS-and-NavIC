function prof = importCustomTrajectory(csvFile, name, opt)
%IMPORTCUSTOMTRAJECTORY Prepare a user trajectory CSV for IF generation.
%   prof = importCustomTrajectory(csvFile, name, opt)
%
%   csvFile - trajectory CSV, one header line, columns
%               Time [sec], Latitude [deg], Longitude [deg], Height [meter],
%               Roll Angle [deg], Pitch Angle [deg], Heading [deg]
%             (any sample rate; attitude = aerospace Euler angles w.r.t. NED)
%   name    - profile name (outputs Fighter_<name>.*)
%   opt     - optional struct: outDir (default <simulator>\FighterProfiles),
%             accelRangeG (default 20), seed (IMU errors, default 1)
%
%   Writes, like generateFighterProfile.m,
%     Fighter_<name>.csv                 1 kHz cubic-spline resampling (ECEF)
%     Fighter_<name>_IMU_200Hz.mat/.csv  iNAT-FSLG-01 IMU data with truth
%   The IMU truth comes from exact derivatives of the same ECEF spline that
%   feeds the IF simulator. The angular rate is formed from Euler-angle
%   rates, which suits smooth attitude histories; aerobatic trajectories
%   (loops, rolls through +/-90 deg pitch) should be built with the custom
%   manoeuvre builder instead, which carries the attitude as a DCM.

simDir = fileparts(mfilename('fullpath'));
if nargin < 3, opt = struct(); end
if ~isfield(opt, 'outDir'),      opt.outDir = fullfile(simDir, 'FighterProfiles'); end
if ~isfield(opt, 'accelRangeG'), opt.accelRangeG = 20; end
if ~isfield(opt, 'seed'),        opt.seed = 1; end
if ~exist(opt.outDir, 'dir'), mkdir(opt.outDir); end

% ---- checks
raw = readmatrix(csvFile, 'NumHeaderLines', 1);
if size(raw, 2) < 7, error('importCustomTrajectory: %s must have 7 columns (time, lat, lon, h, roll, pitch, heading).', csvFile); end
raw = raw(:, 1:7);
if any(~isfinite(raw(:))), error('importCustomTrajectory: %s contains non-numeric or missing values.', csvFile); end
if any(diff(raw(:,1)) <= 0), error('importCustomTrajectory: time column must be strictly increasing.'); end
T = raw(end,1) - raw(1,1);
if T < 13, error('importCustomTrajectory: trajectory spans %.1f s; at least 13 s is needed (receiver start-up).', T); end
if any(abs(raw(:,2)) > 90) || any(raw(:,4) < -500), error('importCustomTrajectory: latitude or height out of range.'); end
if T > 60, warning('importCustomTrajectory: %.0f s trajectory; generateIFdata writes at most 60 s unless opt.ms is set.', T); end

% ---- 1 kHz trajectory and IMU (the simulator expects the trajectory to start at t = 0)
base = fullfile(opt.outDir, ['Fighter_' name]);
tmp  = [tempname '.csv'];  c = onCleanup(@() deleteIf(tmp));
raw(:,1) = raw(:,1) - raw(1,1);
fid = fopen(tmp, 'w');
fprintf(fid, 'Time [sec],Latitude [deg],Longitude [deg],Height [meter],Roll Angle [deg],Pitch Angle [deg],Heading [deg]\n');
fprintf(fid, '%.6f,%.12f,%.12f,%.6f,%.10f,%.10f,%.10f\n', raw');
fclose(fid);
resampleTrajectorySpline(tmp, 1000, [base '.csv']);
imu = generateIMU_iNAT_FSLG01([base '.csv'], 200, opt.accelRangeG, opt.seed, [base '_IMU_200Hz']);

% ---- body-to-ECEF attitude truth (as in generateFighterProfile.m) and profile record
d2r = pi/180;  n = numel(imu.truth.t);  Cbe = zeros(3,3,n);
for k = 1:n
    la = imu.truth.llh(k,1)*d2r;  lo = imu.truth.llh(k,2)*d2r;  e = imu.truth.attRad(k,:);
    Cbe(:,:,k) = cne(la, lo)*cbn(e(1), e(2), e(3));
end
imu.truth.Cbe = Cbe;

% profile statistics in the form written by generateFighterProfile.m (used by
% the receiver summaries); computed from the 200 Hz truth
g0 = 9.80665;  dt = 1/200;
fT = vecnorm(imu.fTrue, 2, 2)/g0;
nz = -imu.fTrue(:,3)/g0;
win = 2;                                                    % 10 ms smoothing
dn  = movmean(gradient(movmean(nz, win), dt), win);
aE  = gradient(imu.truth.vEcef', dt)';
[~, daE] = gradient(movmean(aE, win), 1, dt);
V   = vecnorm(imu.truth.vEcef, 2, 2);
st = struct('model', 'trajectory', 'maxG', max(fT), 'minNz', min(nz), 'maxNz', max(nz), ...
    'maxOnsetGps', max(abs(dn)), 'maxJerkGps', max(vecnorm(movmean(daE, win), 2, 2))/g0, ...
    'maxRollRate', max(abs(imu.wTrue(:,1)))/d2r, 'maxRate', max(vecnorm(imu.wTrue, 2, 2))/d2r, ...
    'hRange', [min(imu.truth.llh(:,3)) max(imu.truth.llh(:,3))], 'VRange', [min(V) max(V)], ...
    'invertedFrac', 0, 'maneuverWin', [min(12, T) T], 'alphaRange', [0 0], 'maxBeta', 0);
inv = false(n,1);                                           % antenna (body -z) below the horizon
for k = 1:n
    Cbn = cbn(imu.truth.attRad(k,1), imu.truth.attRad(k,2), imu.truth.attRad(k,3));
    inv(k) = Cbn(3,3) < 0;
end
st.invertedFrac = mean(inv);
imu.profile = struct('name', name, 'description', ['Custom trajectory from ' csvFile], ...
                     'maneuverWin', st.maneuverWin, 'stats', st);
imu.profile.flight = struct('t', imu.truth.t(:), 'n', nz, 'h', imu.truth.llh(:,3), 'V', V, ...
                            'eulDeg', imu.truth.attRad/d2r, 'wDeg', imu.wTrue/d2r);
save([base '_IMU_200Hz.mat'], 'imu');
prof = imu.profile;
fprintf('Custom trajectory %s: %.1f s, max |f| %.1f g -> %s.csv\n', name, T, max(fT), base);
end

function C = cbn(ph, th, ps)
% body -> NED for roll ph, pitch th, heading ps (3-2-1)
cph = cos(ph); sph = sin(ph); cth = cos(th); sth = sin(th); cps = cos(ps); sps = sin(ps);
C = [cth*cps, sph*sth*cps - cph*sps, cph*sth*cps + sph*sps;
     cth*sps, sph*sth*sps + cph*cps, cph*sth*sps - sph*cps;
     -sth,    sph*cth,               cph*cth];
end

function C = cne(la, lo)
% NED -> ECEF at geodetic latitude la, longitude lo
C = [-sin(la)*cos(lo), -sin(lo), -cos(la)*cos(lo);
     -sin(la)*sin(lo),  cos(lo), -cos(la)*sin(lo);
      cos(la),          0,       -sin(la)];
end

function deleteIf(f)
if exist(f, 'file'), delete(f); end
end
