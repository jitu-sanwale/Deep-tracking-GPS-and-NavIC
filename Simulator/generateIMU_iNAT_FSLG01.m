function imu = generateIMU_iNAT_FSLG01(trajFile, fsIMU, accelRangeG, seed, outBase)
%GENERATEIMU_INAT_FSLG01 Synthesise iMAR iNAT-FSLG-01 IMU data along a trajectory.
%   imu = generateIMU_iNAT_FSLG01(trajFile, fsIMU, accelRangeG, seed, outBase)
%
%   trajFile    - trajectory CSV (signalGenerating.m format). The truth is the
%                 ECEF cubic spline of trajectorySplineECEF.m, i.e. the same
%                 one resampleTrajectorySpline.m feeds to the IF simulator.
%   fsIMU       - output rate [Hz] (default 200, i.e. the unit's 1 kHz
%                 internal rate decimated by 5)
%   accelRangeG - accelerometer range [g]; [] -> prompt (AccelRangeSetting.m,
%                 very high default, 20 g datasheet option)
%   seed        - rng seed for the error realisation (default 1)
%   outBase     - output file base name (default 'iNAT_FSLG01_IMU_<fs>Hz');
%                 writes <outBase>.mat and <outBase>.csv
%
%   Outputs are interval averages over (t-dt, t]:  f = dv/dt [m/s^2] and
%   w = dtheta/dt [rad/s] in the body frame (x fwd, y right, z down, body
%   attitude from the trajectory roll/pitch/heading; all zero for the
%   projectile CSV, so body = local NED there).
%
%   Truth (internal 4 kHz grid, integrated with Simpson's rule to dv,dtheta):
%       f_e    = a_e + 2*Omega_ie*v_e - g_e           (g_e: WGS-84 J2 gravity
%                                                      incl. centrifugal)
%       f_b    = C_n^b C_e^n f_e
%       w_ib^b = C_n^b (w_ie^n + w_en^n) + w_nb^b(Euler rates)
%   The Euler-rate term assumes smooth attitude columns; for aerobatic
%   profiles use generateFighterProfile.m, which carries the attitude as a
%   DCM. Sensor errors: imuErrorModel_iNAT.m (full datasheet model).

if nargin < 2 || isempty(fsIMU), fsIMU = 200; end
if nargin < 3, accelRangeG = []; end
if nargin < 4 || isempty(seed), seed = 1; end
if nargin < 5 || isempty(outBase), outBase = sprintf('iNAT_FSLG01_IMU_%dHz', fsIMU); end

accelRangeG = AccelRangeSetting(accelRangeG);

%% Constants (WGS-84)
g0    = 9.80665;
omgIE = 7.292115e-5;
R0    = 6378137;  ecc = 0.0818191908425;
d2r   = pi/180;

%% Truth on the internal grid
traj  = trajectorySplineECEF(trajFile);
dtIMU = 1/fsIMU;
nSub  = 20;                              % Simpson sub-intervals per output (4 kHz)
nOut  = floor(round((traj.t1 - traj.t0)/dtIMU, 6));
tF    = traj.t0 + (0:nOut*nSub)'*(dtIMU/nSub);
tF(end) = min(tF(end), traj.t1);

r = ppval(traj.ppPos, tF')';  v = ppval(traj.ppVel, tF')';  a = ppval(traj.ppAcc, tF')';
att  = ppval(traj.ppAtt, tF')';  attR = ppval(traj.ppAttRate, tF')';
[latD, ~, h] = ecef2geodeticWGS84(r);
lat = latD*d2r;  lon = atan2(r(:,2), r(:,1));

nF = numel(tF);
fB = zeros(nF,3);  wB = zeros(nF,3);
OmgIE = [0 -omgIE 0; omgIE 0 0; 0 0 0];
for k = 1:nF
    fE = a(k,:)' + 2*OmgIE*v(k,:)' - gravityECEFsim(r(k,:)');

    sL = sin(lat(k)); cL = cos(lat(k)); sl = sin(lon(k)); cl = cos(lon(k));
    Cen = [-sL*cl -sL*sl cL; -sl cl 0; -cL*cl -cL*sl -sL];
    vN  = Cen*v(k,:)';
    RN  = R0*(1-ecc^2)/(1-ecc^2*sL^2)^1.5;  RE = R0/sqrt(1-ecc^2*sL^2);
    wIEn = omgIE*[cL; 0; -sL];
    wENn = [vN(2)/(RE+h(k)); -vN(1)/(RN+h(k)); -vN(2)*tan(lat(k))/(RE+h(k))];

    ph = att(k,1); th = att(k,2); ps = att(k,3);
    Cnb = eulerToCnb(ph, th, ps);
    dph = attR(k,1); dth = attR(k,2); dps = attR(k,3);
    wNBb = [dph - dps*sin(th); dth*cos(ph) + dps*sin(ph)*cos(th); -dth*sin(ph) + dps*cos(ph)*cos(th)];

    fB(k,:) = (Cnb*Cen*fE)';
    wB(k,:) = (Cnb*(wIEn + wENn) + wNBb)';
end

%% Sensor errors
imu = imuErrorModel_iNAT(tF, fB, wB, fsIMU, nSub, accelRangeG, seed);

%% Truth at the output instants
tq  = imu.t';
rT  = ppval(traj.ppPos, tq)';  vT = ppval(traj.ppVel, tq)';
attT = ppval(traj.ppAtt, tq)';
[latT, lonT, hT] = ecef2geodeticWGS84(rT);
imu.truth.t   = imu.t;  imu.truth.rEcef = rT;  imu.truth.vEcef = vT;
imu.truth.llh = [latT lonT hT];  imu.truth.attRad = attT;
imu.trajFile  = trajFile;

save([outBase '.mat'], 'imu');
writeImuCsv([outBase '.csv'], imu);

fMag = vecnorm(imu.fTrue,2,2)/g0;
fprintf('iNAT-FSLG-01 IMU: %d samples @ %g Hz, %.3f-%.3f s -> %s.mat/.csv\n', numel(imu.t), fsIMU, imu.t(1), imu.t(end), outBase);
fprintf('  true |f| max %.1f g; accel range +/-%g g; %d output samples contain clipping\n', ...
    max(fMag), accelRangeG, nnz(imu.clipped));
fprintf('  turn-on bias: gyro %s deg/h, accel %s ug\n', ...
    mat2str(round(imu.err.bg0'/d2r*3600,4)), mat2str(round(imu.err.ba0'/g0*1e6,1)));
fprintf('  rms(meas-true): f %s m/s^2, w %s deg/h\n', ...
    mat2str(round(rms(imu.f-imu.fTrue),5)), mat2str(round(rms(imu.w-imu.wTrue)/d2r*3600,4)));
end

function C = eulerToCnb(ph, th, ps)
% Groves Euler_to_CTM: coordinate transformation matrix n -> b
sp = sin(ph); cp = cos(ph); st = sin(th); ct = cos(th); ss = sin(ps); cs = cos(ps);
C = [ct*cs,               ct*ss,               -st;
     -cp*ss + sp*st*cs,   cp*cs + sp*st*ss,    sp*ct;
     sp*ss + cp*st*cs,   -sp*cs + cp*st*ss,    cp*ct];
end
