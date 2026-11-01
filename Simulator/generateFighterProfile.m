function prof = generateFighterProfile(name, opt)
%GENERATEFIGHTERPROFILE Fighter-level manoeuvre trajectory with attitude + IMU.
%   prof = generateFighterProfile(name, opt)
%
%   name - one of fighterProfileList() (e.g. 'turn', 'break_gonset', ...)
%   opt  - optional struct: outDir (default 'FighterProfiles'), accelRangeG
%          (default 20 = iNAT-FSLG-01 datasheet; fighter loads stay below it),
%          seed (IMU error realisation, default 1), T (duration, default 60 s)
%
%   Each profile (fighterPhases.m) is a sequence of phases giving a roll-rate
%   reference p, a normal-load-factor reference n and dV/dt. Two flight
%   models fly them (P.model):
%     '6dof'      - airborne profiles: nonlinear 6-DOF F-16 (NASA TP-1538
%                   aerodynamics, MATLAB port of F16Model.jl) under a
%                   dynamic-inversion autopilot (F16Aircraft/integrate6dof.m): real angle
%                   of attack, sideslip, roll/pitch dynamics, thrust and drag.
%     'pointmass' - ground-contact profiles (catapult, landing), which the
%                   6-DOF model cannot represent (no landing gear): 3-DOF point
%                   mass, coordinated flight, q = (n g - g_bz)/V, r = g_by/V,
%                   g-onset limited to P.jerkGps, roll acceleration to
%                   P.pddMax; ground phases freeze the attitude.
%   Both integrate the attitude as a DCM (RK4, 4 kHz), so loops, rolls and
%   inverted flight have no Euler singularity.
%
%   Outputs (in opt.outDir):
%     Fighter_<name>.csv          1 kHz trajectory for signalGenerating.m
%                                 (lat, lon, h, roll, pitch, heading; standard
%                                 aerospace Euler angles w.r.t. local NED)
%     Fighter_<name>_IMU_200Hz.mat/.csv  iNAT-FSLG-01 data (imuErrorModel_iNAT)
%                                 with truth (ECEF position/velocity, C_b^e)
%   The IMU truth is computed exactly in ECEF from the integrated path:
%       f_b = C_e^b (a_e + 2 W_ie v_e - g_e),  w_ib^b = C_e^b w_ie + w_nb^b
%   (the tangent frame is fixed in ECEF), so it is consistent with the IF
%   signal to the integration accuracy. P.vibe(t), if given, adds airframe
%   vibration to the body specific force (no measurable CG displacement).

if nargin < 2, opt = struct(); end
if ~isfield(opt,'outDir'),      opt.outDir = 'FighterProfiles'; end
if ~isfield(opt,'accelRangeG'), opt.accelRangeG = 20; end
if ~isfield(opt,'seed'),        opt.seed = 1; end
if ~isfield(opt,'T'),           opt.T = 60; end
if ~exist(opt.outDir,'dir'), mkdir(opt.outDir); end

g0 = 9.80665;  d2r = pi/180;
fsIMU = 200;  nSub = 20;  dt = 1/(fsIMU*nSub);         % 4 kHz internal grid

%% Profile definition and flight
P = fighterPhases(name, opt);
lat0 = 20.1085;  lon0 = 73.9532;
if strcmp(P.model, '6dof')
    addpath(fullfile(fileparts(mfilename('fullpath')), 'F16Aircraft'));
    % 6-DOF integrated at 1 kHz (well above its short-period/roll-mode and
    % 20 ms actuator dynamics; the 1 kHz nodes are exactly the CSV the IF
    % simulator uses), then spline-resampled to the 4 kHz IMU truth grid
    [rec, phStart, iph] = integrate6dof(P, opt.T, 1e-3, g0);
    rec = resampleRec(rec, dt);
else
    [rec, phStart, iph] = integratePointMass(P, opt.T, dt, g0);
end
if iph < numel(P.phases)
    warning('generateFighterProfile:%s ended in phase %d of %d (%s) - the manoeuvre did not complete in %g s.', ...
        name, iph, numel(P.phases), P.phases(iph).name, opt.T);
end
N = numel(rec.t);

%% ECEF truth
[x0, y0, z0] = geo2cartd(lat0, lon0, 0, 5);
re0 = [x0; y0; z0];
Cne = cen(lat0*d2r, lon0*d2r)';
omgIE = 7.292115e-5;  OmgIE = [0 -omgIE 0; omgIE 0 0; 0 0 0];
rE = (re0 + Cne*rec.r')';  vE = (Cne*rec.v')';  aE = (Cne*rec.a')';
fB = zeros(N,3);  wB = zeros(N,3);  Cbe = zeros(3,3,N);
for k = 1:N
    Cbe(:,:,k) = Cne*rec.C(:,:,k);
    fe = aE(k,:)' + 2*OmgIE*vE(k,:)' - gravityECEFsim(rE(k,:)');
    fB(k,:) = (Cbe(:,:,k)'*fe)';
    wB(k,:) = (Cbe(:,:,k)'*[0;0;omgIE] + rec.w(k,:)')';
end
if isfield(P,'vibe') && ~isempty(P.vibe)
    fB = fB + P.vibe(rec.t);
end
[latD, lonD, hE] = ecef2geodeticWGS84(rE);
eul = zeros(N,3);
for k = 1:N
    Cbn = cen(latD(k)*d2r, lonD(k)*d2r)*Cbe(:,:,k);       % body -> local NED
    eul(k,:) = [atan2(Cbn(3,2), Cbn(3,3)), -asin(max(min(Cbn(3,1),1),-1)), atan2(Cbn(2,1), Cbn(1,1))];
end

%% 1 kHz trajectory CSV
i1k = 1:4:N;
base = fullfile(opt.outDir, ['Fighter_' name]);
fid = fopen([base '.csv'], 'w');
fprintf(fid, 'Time [sec],Latitude [deg],Longitude [deg],Height [meter],Roll Angle [deg],Pitch Angle [deg],Heading [deg]\n');
fprintf(fid, '%.6f,%.12f,%.12f,%.6f,%.6f,%.6f,%.6f\n', ...
    [rec.t(i1k) latD(i1k) lonD(i1k) hE(i1k) eul(i1k,1)/d2r eul(i1k,2)/d2r mod(eul(i1k,3)/d2r,360)]');
fclose(fid);

%% IMU
imu = imuErrorModel_iNAT(rec.t, fB, wB, fsIMU, nSub, opt.accelRangeG, opt.seed);
iI = nSub+1:nSub:N;
imu.truth.t      = imu.t;
imu.truth.rEcef  = rE(iI,:);   imu.truth.vEcef = vE(iI,:);
imu.truth.llh    = [latD(iI) lonD(iI) hE(iI)];
imu.truth.attRad = eul(iI,:);
imu.truth.Cbe    = Cbe(:,:,iI);
imu.trajFile     = [base '.csv'];

%% profile statistics (from the truth)
nz   = -fB(:,3)/g0;                                     % normal load factor
win  = round(0.01/dt);                                   % 10 ms smoothing for rates of change
dn   = movmean(gradient(movmean(nz, win), dt), win);
[~, daE] = gradient(movmean(aE, win), 1, dt);           % d/dt along rows (time)
jerk = vecnorm(movmean(daE, win), 2, 2)/g0;
iFirst = find([P.phases.maneuver], 1);  iLast = find([P.phases.maneuver], 1, 'last');
mw = [phStart(iFirst), opt.T];                           % manoeuvre window [s]
if iLast < numel(P.phases) && phStart(iLast+1) > mw(1), mw(2) = phStart(iLast+1); end
st.model        = P.model;
st.maxG         = max(vecnorm(fB,2,2))/g0;
st.minNz        = min(nz);  st.maxNz = max(nz);
st.maxOnsetGps  = max(abs(dn));
st.maxJerkGps   = max(jerk);
st.maxRollRate  = max(abs(wB(:,1)))/d2r;
st.maxRate      = max(vecnorm(wB,2,2))/d2r;
st.hRange       = [min(hE) max(hE)];
st.VRange       = [min(vecnorm(rec.v,2,2)) max(vecnorm(rec.v,2,2))];
st.invertedFrac = mean(squeeze(rec.C(3,3,:)) < 0);      % antenna (body -z) facing below the horizon
st.maneuverWin  = mw;
st.alphaRange   = [0 0];  st.maxBeta = 0;
if isfield(rec,'alpha'), st.alphaRange = [min(rec.alpha) max(rec.alpha)]; st.maxBeta = max(abs(rec.beta)); end

imu.profile.name        = name;
imu.profile.description = P.description;
imu.profile.maneuverWin = mw;
imu.profile.stats       = st;
imu.profile.phaseNames  = {P.phases.name};
imu.profile.phaseStart  = phStart;
i100 = 1:40:N;                                           % 100 Hz flight record
imu.profile.flight.t = rec.t(i100);  imu.profile.flight.n = nz(i100);
imu.profile.flight.h = hE(i100);     imu.profile.flight.V = vecnorm(rec.v(i100,:),2,2);
imu.profile.flight.eulDeg = eul(i100,:)/d2r;  imu.profile.flight.wDeg = wB(i100,:)/d2r;
if isfield(rec,'alpha')
    imu.profile.flight.alpha = rec.alpha(i100);  imu.profile.flight.beta = rec.beta(i100);
    imu.profile.flight.ctrl  = rec.ctrl(i100,:);
end
save([base '_IMU_200Hz.mat'], 'imu');
writeImuCsv([base '_IMU_200Hz.csv'], imu);

prof = imu.profile;
fprintf(['%-17s [%s] %s\n   max %.1f g (n %.1f..%.1f), onset %.1f g/s, jerk %.1f g/s, roll %.0f deg/s, |w| %.0f deg/s, ' ...
         'alpha %.1f..%.1f deg, |beta| %.1f deg\n   h %.0f-%.0f m, V %.0f-%.0f m/s, inverted %.0f %%, manoeuvre %.1f-%.1f s, clipped %d\n'], ...
    name, P.model, P.description, st.maxG, st.minNz, st.maxNz, st.maxOnsetGps, st.maxJerkGps, st.maxRollRate, st.maxRate, ...
    st.alphaRange, st.maxBeta, st.hRange, st.VRange, 100*st.invertedFrac, mw, nnz(imu.clipped));
end

%% ========================================================================
function r4 = resampleRec(rec, dt)
% spline-resample a 1 kHz flight record onto the dt grid (DCM: linear
% interpolation + re-orthonormalisation; < 0.3 deg rotation per 1 ms step)
t1 = rec.t;  t4 = (0:round(t1(end)/dt))'*dt;
r4.t = t4;
for f = {'r','v','a','w','n','alpha','beta','ctrl','V'}
    r4.(f{1}) = interp1(t1, rec.(f{1}), t4, 'spline');
end
r4.ph = interp1(t1, rec.ph, t4, 'previous');
C9 = interp1(t1, reshape(rec.C, 9, [])', t4, 'linear');
r4.C = zeros(3,3,numel(t4));
for k = 1:numel(t4)
    C = reshape(C9(k,:), 3, 3);
    r4.C(:,:,k) = C*(3*eye(3) - C'*C)/2;
end
end

function [rec, phStart, iph] = integratePointMass(P, T, dt, g0)
d2r = pi/180;
N  = round(T/dt) + 1;
x  = zeros(17,1);
C0 = cnb(0, 0, P.hdg0*d2r)';
x(1:3) = [0; 0; -P.h0];  x(4) = P.V0;  x(5:13) = C0(:);  x(16) = 1;
rec.t = (0:N-1)'*dt;
rec.r = zeros(N,3); rec.v = zeros(N,3); rec.a = zeros(N,3);
rec.C = zeros(3,3,N); rec.w = zeros(N,3); rec.n = zeros(N,1); rec.ph = zeros(N,1);
iph = 1;  tPh = 0;  s0 = stateOf(x);  phStart = zeros(1,numel(P.phases));
for k = 1:N
    t = rec.t(k);
    [~, ex] = deriv(x, P.phases(iph), t - tPh, s0, P, g0);
    C = reshape(x(5:13),3,3);
    rec.r(k,:) = x(1:3)';  rec.v(k,:) = (C(:,1)*x(4))';  rec.a(k,:) = (C*ex.ab)';
    rec.C(:,:,k) = C;  rec.w(k,:) = ex.w';  rec.n(k) = ex.n;  rec.ph(k) = iph;
    if k == N, break; end
    ph = P.phases(iph);  tau = t - tPh;
    k1 = deriv(x,           ph, tau,        s0, P, g0);
    k2 = deriv(x + dt/2*k1, ph, tau + dt/2, s0, P, g0);
    k3 = deriv(x + dt/2*k2, ph, tau + dt/2, s0, P, g0);
    k4 = deriv(x + dt*k3,   ph, tau + dt,   s0, P, g0);
    x  = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    C  = reshape(x(5:13),3,3);  [U,~,V] = svd(C);  C = U*V';  x(5:13) = C(:);   % re-orthonormalise
    if x(4) < 0, x(4) = 0; end
    s = stateOf(x);
    if iph < numel(P.phases) && P.phases(iph).done(t + dt - tPh, s, s0)
        iph = iph + 1;  tPh = t + dt;  s0 = s;  x(14:15) = 0;  phStart(iph) = tPh;
    end
end
end

function [dx, ex] = deriv(x, ph, tau, s0, P, g0)
s = stateOf(x);
[pRef, nRef, Vdot, ground, nDir] = ph.cmd(tau, s, s0);
nAct = x(16);  pAct = x(17);
dn = max(min((nRef - nAct)/0.1, P.jerkGps), -P.jerkGps);
dp = max(min((pRef - pAct)/0.05, P.pddMax), -P.pddMax);
C = s.C;  V = s.V;
if ground
    p = 0; q = 0; r = 0; n = 1;
    if V <= 0 && Vdot < 0, Vdot = 0; end
else
    n = nAct + nDir;  p = pAct;
    Ve = max(V, 30);
    gb = g0*C(3,:)';
    q = (n*g0 - gb(3))/Ve;
    r = gb(2)/Ve;
end
w  = [p; q; r];
dC = C*skew(w);
dx = [C(:,1)*V; Vdot; dC(:); q; p; dn; dp];
ex.w = w;  ex.n = n;  ex.ab = [Vdot; r*V; -q*V];
end

function s = stateOf(x)
s.r = x(1:3);  s.V = x(4);  s.C = reshape(x(5:13),3,3);
s.qint = x(14);  s.pint = x(15);  s.nAct = x(16);  s.alpha = 0;  s.w = [0; 0; 0];
C = s.C;
s.phi   = atan2(C(3,2), C(3,3));
s.theta = -asin(max(min(C(3,1),1),-1));
s.gamma = s.theta;                          % zero angle of attack
s.psi   = atan2(C(2,1), C(1,1));
s.h     = -x(3);
s.upright = C(3,3) > 0;
end

function C = cnb(ph, th, ps)
sp = sin(ph); cp = cos(ph); st = sin(th); ct = cos(th); ss = sin(ps); cs = cos(ps);
C = [ct*cs, ct*ss, -st; -cp*ss+sp*st*cs, cp*cs+sp*st*ss, sp*ct; sp*ss+cp*st*cs, -sp*cs+cp*st*ss, cp*ct];
end

function C = cen(lat, lon)
sL = sin(lat); cL = cos(lat); sl = sin(lon); cl = cos(lon);
C = [-sL*cl -sL*sl cL; -sl cl 0; -cL*cl -cL*sl -sL];
end

function S = skew(a)
S = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
end
