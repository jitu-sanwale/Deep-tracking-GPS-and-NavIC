function traj = trajectorySplineECEF(trajFile)
%TRAJECTORYSPLINEECEF Cubic-spline model of a trajectory CSV, built in ECEF.
%   traj = trajectorySplineECEF(trajFile)
%
%   The CSV (Time [sec], Latitude [deg], Longitude [deg], Height [meter],
%   Roll/Pitch/Heading [deg]) is converted to WGS-84 ECEF and one
%   not-a-knot cubic spline is fitted per ECEF axis. This single spline is
%   the common ground truth for
%     - the 1 kHz resampled trajectory CSV fed to signalGenerating.m
%       (resampleTrajectorySpline.m), and
%     - the IMU truth (generateIMU_iNAT_FSLG01.m), which takes velocity and
%       acceleration as exact derivatives of the same spline, so there is no
%       numerical differentiation noise (a second difference of ECEF
%       coordinates at 1 ms would carry ~1e-3 m/s^2 of round-off).
%   The spline is fitted in ECEF rather than in lat/lon/h so that its
%   derivatives are directly the ECEF velocity/acceleration.
%
%   Output struct fields:
%       t0, t1          - time span [s]
%       ppPos           - 3-vector pp form, r_e(t) [m]
%       ppVel, ppAcc    - its first and second derivatives
%       ppAtt           - 3-vector pp form, [roll pitch heading](t) [rad]
%                         (heading unwrapped before fitting)
%       ppAttRate       - derivative of ppAtt [rad/s]
%       raw             - the CSV contents

raw = readmatrix(trajFile, 'NumHeaderLines', 1);
t   = raw(:,1);
if size(raw,2) >= 7
    att = raw(:,5:7) * pi/180;
else
    att = zeros(numel(t), 3);
end
att(:,3) = unwrap(att(:,3));

xyz = zeros(numel(t), 3);
for k = 1:numel(t)
    [xyz(k,1), xyz(k,2), xyz(k,3)] = geo2cartd(raw(k,2), raw(k,3), raw(k,4), 5); % WGS-84
end

traj.t0        = t(1);
traj.t1        = t(end);
traj.ppPos     = spline(t', xyz');
traj.ppVel     = fnder(traj.ppPos, 1);
traj.ppAcc     = fnder(traj.ppPos, 2);
traj.ppAtt     = spline(t', att');
traj.ppAttRate = fnder(traj.ppAtt, 1);
traj.raw       = raw;
end
