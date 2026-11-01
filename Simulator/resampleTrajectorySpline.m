function outFile = resampleTrajectorySpline(inFile, fsOut, outFile)
%RESAMPLETRAJECTORYSPLINE Resample a trajectory CSV with a cubic spline.
%   outFile = resampleTrajectorySpline(inFile, fsOut, outFile)
%
%   inFile  - trajectory CSV (signalGenerating.m format)
%   fsOut   - output rate [Hz] (default 1000, the simulator's 1 ms block
%             grid, so signalGenerating.m's own linear interpolation becomes
%             an identity and the IF signal follows the spline exactly)
%   outFile - output CSV (default <inFile>_<fsOut/1000>kHz.csv)
%
%   The spline is fitted in ECEF (trajectorySplineECEF.m) and evaluated on
%   the new grid, then converted back to WGS-84 lat/lon/h. Attitude columns
%   are resampled with their own spline.

if nargin < 2 || isempty(fsOut), fsOut = 1000; end
if nargin < 3 || isempty(outFile)
    [~, base] = fileparts(inFile);
    outFile = sprintf('%s_%gkHz.csv', base, fsOut/1000);
end

traj = trajectorySplineECEF(inFile);
nOut = floor(round((traj.t1 - traj.t0)*fsOut, 6)) + 1;
tq   = traj.t0 + (0:nOut-1)'/fsOut;

xyz = ppval(traj.ppPos, tq')';
att = ppval(traj.ppAtt, tq')' * 180/pi;
att(:,3) = mod(att(:,3), 360);
att(abs(att) < 1e-12) = 0;
[lat, lon, h] = ecef2geodeticWGS84(xyz);

fid = fopen(outFile, 'w');
fprintf(fid, 'Time [sec],Latitude [deg],Longitude [deg],Height [meter],Roll Angle [deg],Pitch Angle [deg],Heading [deg]\n');
fprintf(fid, '%.6f,%.12f,%.12f,%.6f,%.10f,%.10f,%.10f\n', [tq lat lon h att]');
fclose(fid);

% Report how far the spline departs from the original samples' linear
% interpolation (i.e. what the simulator would have used before).
raw  = traj.raw;
xyzL = zeros(nOut,3);
xyzR = zeros(size(raw,1),3);
for k = 1:size(raw,1)
    [xyzR(k,1), xyzR(k,2), xyzR(k,3)] = geo2cartd(raw(k,2), raw(k,3), raw(k,4), 5);
end
for ax = 1:3, xyzL(:,ax) = interp1(raw(:,1), xyzR(:,ax), tq, 'linear'); end
dev = vecnorm(xyz - xyzL, 2, 2);
fprintf('Resampled %s: %d -> %d points (%.0f Hz), %.3f-%.3f s -> %s\n', ...
    inFile, size(raw,1), nOut, fsOut, tq(1), tq(end), outFile);
fprintf('  spline vs. linear interpolation: max %.3f m, rms %.3f m\n', max(dev), sqrt(mean(dev.^2)));
end
