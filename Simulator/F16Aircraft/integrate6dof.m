function [rec, phStart, iph] = integrate6dof(P, T, dt, g0)
%INTEGRATE6DOF Fly a manoeuvre profile with a nonlinear 6-DOF aircraft model.
%   [rec, phStart, iphEnd] = integrate6dof(P, T, dt, g0)
%
%   Airframe (P.aircraft, see aircraftConfig below):
%     'F16' - NASA TP-1538 / Stevens & Lewis F-16 (F16Aircraft\f16Aero.m, a
%             MATLAB port of F16Model.jl), default
%   Rigid-body equations as in Stengel's FLIGHT v2 EoMQv2.m / Stevens & Lewis
%   (flat-earth NED, body-axis velocity, Euler moment equations with Ixz);
%   attitude carried as a DCM (no Euler singularity in loops). RK4 at dt.
%
%   Flight control (1 kHz, zero-order hold) - nonlinear dynamic inversion
%   on the aircraft's own aero model, tracking the phase references of
%   fighterPhases.m:
%     pitch : load factor n_ref -> q_cmd = (n g - g_bz + p v)/u + PI on n error
%     roll  : roll rate p_ref
%     yaw   : sideslip -> 0 (coordinated)
%     speed : dV/dt reference -> thrust
%   Desired angular accelerations are converted to elevator / aileron /
%   rudder through the moment equations and the model's control
%   effectiveness. Surface limits and rates per aircraft, first-order 20 ms
%   actuators, thrust lag. A phase may override any control (6th output of
%   its command function, NaN = keep the autopilot) - used for the
%   departure and the snap roll. P.wind(t) (NED gust velocity, m/s) and
%   P.dCl(t) (rolling-moment coefficient disturbance, e.g. a wake vortex)
%   are optional.
%
%   rec: t, r (NED pos), v (NED vel), a (NED accel), C (C_b^n), w (body
%   rates), n (normal load factor), ph (phase index), alpha, beta (deg),
%   ctrl (elevator, aileron, rudder deg, thrust fraction), V (airspeed)

if ~isfield(P,'aircraft') || isempty(P.aircraft), P.aircraft = 'F16'; end
ac = aircraftConfig(P.aircraft);
Kq = 12; Kp = 8; Kr = 10; Kb = 5; Kn = 1.0; Kni = 1.5;   % roll gain matched to the slow (80 deg/s) aileron
nCtl = max(1, round(1e-3/dt));                     % control at 1 kHz
if ~isfield(P,'wind') || isempty(P.wind), P.wind = @(t) [0;0;0]; end
if ~isfield(P,'dCl')  || isempty(P.dCl),  P.dCl  = @(t) 0; end

%% trim (straight and level at P.V0, P.h0)
[ab0, de0, th0] = trimLevel(ac, P.V0, P.h0, g0);
C0 = cnbE(0, ab0, P.hdg0*pi/180)';
x = zeros(24,1);
x(1:3) = [0; 0; -P.h0];  x(4:6) = P.V0*[cos(ab0); 0; sin(ab0)];
x(7:15) = C0(:);  x(21:24) = [de0; 0; 0; th0];
dcmd = x(21:24);  nInt = 0;
fprintf('  6-DOF %s trim: alpha %.2f deg, elevator %.2f deg, thrust %.3f of max\n', ac.name, ab0*180/pi, de0*180/pi, th0);

N = round(T/dt) + 1;
rec.t = (0:N-1)'*dt;
rec.r = zeros(N,3); rec.v = zeros(N,3); rec.a = zeros(N,3); rec.C = zeros(3,3,N);
rec.w = zeros(N,3); rec.n = zeros(N,1); rec.ph = zeros(N,1);
rec.alpha = zeros(N,1); rec.beta = zeros(N,1); rec.ctrl = zeros(N,4); rec.V = zeros(N,1);
iph = 1;  tPh = 0;  phStart = zeros(1,numel(P.phases));
[~, ex] = deriv6(x, dcmd, 0, P, g0, ac);
s0 = stateOf6(x, ex);
for k = 1:N
    t = rec.t(k);
    [~, ex] = deriv6(x, dcmd, t, P, g0, ac);
    C = reshape(x(7:15),3,3);
    rec.r(k,:) = x(1:3)';  rec.v(k,:) = (C*x(4:6))';  rec.a(k,:) = (C*ex.ab)';
    rec.C(:,:,k) = C;  rec.w(k,:) = x(16:18)';  rec.n(k) = ex.nz;  rec.ph(k) = iph;
    rec.alpha(k) = ex.alpha*180/pi;  rec.beta(k) = ex.beta*180/pi;  rec.V(k) = ex.V;
    rec.ctrl(k,:) = [x(21:23)'*180/pi, x(24)];
    if k == N, break; end

    % flight control
    if mod(k-1, nCtl) == 0
        s = stateOf6(x, ex);
        [pRef, nRef, Vdot, ~, ~, ovr] = P.phases(iph).cmd(t - tPh, s, s0);
        [dcmd, nInt] = ndi(x, ex, pRef, nRef, Vdot, ovr, nInt, nCtl*dt, g0, ac, ...
                           Kq, Kp, Kr, Kb, Kn, Kni, P.dCl(t), t);
    end

    k1 = deriv6(x,           dcmd, t,        P, g0, ac);
    k2 = deriv6(x + dt/2*k1, dcmd, t + dt/2, P, g0, ac);
    k3 = deriv6(x + dt/2*k2, dcmd, t + dt/2, P, g0, ac);
    k4 = deriv6(x + dt*k3,   dcmd, t + dt,   P, g0, ac);
    x  = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    C  = reshape(x(7:15),3,3);  [U,~,W] = svd(C);  C = U*W';  x(7:15) = C(:);
    if -x(3) < 50
        error('integrate6dof:ground', 'Aircraft below 50 m at t = %.2f s in phase %s.', t, P.phases(iph).name);
    end

    [~, ex] = deriv6(x, dcmd, t + dt, P, g0, ac);
    s = stateOf6(x, ex);
    if iph < numel(P.phases) && P.phases(iph).done(t + dt - tPh, s, s0)
        iph = iph + 1;  tPh = t + dt;  s0 = s;  x(19:20) = 0;  phStart(iph) = tPh;  nInt = 0;
    end
end
end

%% ------------------------------------------------------------------------
function ac = aircraftConfig(name)
d2r = pi/180;
switch upper(name)
    case 'F16'      % F16Model.jl README actuator limits
        ac.name = 'F-16';  ac.aero = @f16Aero;
        addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'F16Aircraft'));
        ac.lim = [25 21.5 30]*d2r;  ac.rate = [60 80 120]*d2r;
        ac.thrMin = 1000/19000;  ac.thrMax = 1;  ac.rateT = 10000/19000;  ac.tauT = 0.2;
    otherwise
        error('integrate6dof: unknown aircraft "%s".', name);
end
ac.tauS = 0.02;
end

function [dx, ex] = deriv6(x, dcmd, t, P, g0, ac, withDeriv)
if nargin < 7, withDeriv = false; end
vb = x(4:6);  C = reshape(x(7:15),3,3);  w = x(16:18);  dl = x(21:24);
h  = -x(3);
Va = vb - C'*P.wind(t);
V  = norm(Va);
al = atan(Va(3)/abs(Va(1)));  be = asin(Va(2)/V);
A  = ac.aero(al, be, V, w, dl, h, withDeriv);
qS = 0.5*A.rho*V^2*A.S;
f  = [(A.CX*qS + A.Thrust)/A.m; A.CY*qS/A.m; A.CZ*qS/A.m];   % specific force (body)
Lb = (A.Cl + P.dCl(t))*qS*A.b;  Mb = A.Cm*qS*A.cBar;  Nb = A.Cn*qS*A.b;
gb = g0*C(3,:)';
dv = f + gb - cross(w, vb);
[p, q, r] = deal(w(1), w(2), w(3));
Ixx = A.Ixx; Iyy = A.Iyy; Izz = A.Izz; Ixz = A.Ixz;  D = Ixx*Izz - Ixz^2;
dp = (Izz*Lb + Ixz*Nb - (Ixz*(Iyy - Ixx - Izz)*p + (Ixz^2 + Izz*(Izz - Iyy))*r)*q)/D;
dq = (Mb - (Ixx - Izz)*p*r - Ixz*(p^2 - r^2))/Iyy;
dr = (Ixz*Lb + Ixx*Nb + (Ixz*(Iyy - Ixx - Izz)*r + (Ixz^2 + Ixx*(Ixx - Iyy))*p)*q)/D;
dC = C*[0 -r q; r 0 -p; -q p 0];
ddl = zeros(4,1);
for i = 1:3
    ddl(i) = max(min((dcmd(i) - dl(i))/ac.tauS, ac.rate(i)), -ac.rate(i));
end
ddl(4) = max(min((dcmd(4) - dl(4))/ac.tauT, ac.rateT), -ac.rateT);
dx = [C*vb; dv; dC(:); dp; dq; dr; q; p; ddl];
ex.f = f;  ex.ab = dv + cross(w, vb);  ex.nz = -f(3)/g0;  ex.alpha = al;  ex.beta = be;
ex.V = V;  ex.Va = Va;  ex.A = A;  ex.qS = qS;  ex.gb = gb;
end

function [dcmd, nInt] = ndi(x, ex, pRef, nRef, Vdot, ovr, nInt, dtc, g0, ac, Kq, Kp, Kr, Kb, Kn, Kni, dClext, t) %#ok<INUSD>
w = x(16:18);  [p, q, r] = deal(w(1), w(2), w(3));  dl = x(21:24);
A  = ac.aero(ex.alpha, ex.beta, ex.V, w, dl, -x(3), true);   % with control effectiveness
Va = ex.Va;  u = max(Va(1), 30);  f = ex.f;  gb = ex.gb;  qS = ex.qS;
% F-16-FLCS-like protections: roll about the velocity vector (stability axis)
% so a roll at high alpha does not convert alpha into sideslip, roll-rate
% limit reduced above 15 deg alpha, alpha limiter at 25 deg / -10 deg
al = ex.alpha;  d2r = pi/180;
pLimFrac = 1 - 0.7*min(max((al - 15*d2r)/(10*d2r), 0), 1);
pRef = max(min(pRef, 6*pLimFrac), -6*pLimFrac);
pStab = pRef;
% load factor -> pitch rate command
nInt = max(min(nInt + (nRef - ex.nz)*dtc, 2), -2);
qc = (nRef*g0 - gb(3))/u + (Kn*(nRef - ex.nz) + Kni*nInt)*g0/u;
if al > 25*d2r,  qc = min(qc, q - 5*(al - 25*d2r));  nInt = min(nInt, 0); end
if al < -10*d2r, qc = max(qc, q - 5*(al + 10*d2r));  nInt = max(nInt, 0); end
rc = (f(2) + gb(2) + p*Va(3) + Kb*Va(2))/u + pStab*sin(al);   % coordination + stability-axis roll
dpD = Kp*(pStab*cos(al) - p);  dqD = Kq*(qc - q);  drD = Kr*(rc - r);
Ixx = A.Ixx; Iyy = A.Iyy; Izz = A.Izz; Ixz = A.Ixz;  D = Ixx*Izz - Ixz^2;
c7 = (Ixz*(Iyy - Ixx - Izz)*p + (Ixz^2 + Izz*(Izz - Iyy))*r)*q;
c9 = (Ixz*(Iyy - Ixx - Izz)*r + (Ixz^2 + Ixx*(Ixx - Iyy))*p)*q;
LN = [Izz Ixz; Ixz Ixx] \ [D*dpD + c7; D*drD - c9];
M  = Iyy*dqD + (Ixx - Izz)*p*r + Ixz*(p^2 - r^2);
dCl = LN(1)/(qS*A.b) - dClext - A.Cl;
dCn = LN(2)/(qS*A.b) - A.Cn;
dCm = M/(qS*A.cBar) - A.Cm;
dcmd = dl;
if abs(A.dCm_de) > 0.02, dcmd(1) = dl(1) + dCm/A.dCm_de; end
B = [A.dCl_da A.dCl_dr; A.dCn_da A.dCn_dr];
if rcond(B) > 1e-6, dar = B \ [dCl; dCn]; dcmd(2:3) = dl(2:3) + dar; end
CX = A.CX;
dcmd(4) = (A.m*(Vdot - gb(1) + q*Va(3) - r*Va(2)) - CX*qS)/(A.StaticThrust*max(A.lapse, 1e-3));
if ~isempty(ovr)
    use = ~isnan(ovr(:));  dcmd(use) = ovr(use);
end
dcmd(1:3) = max(min(dcmd(1:3), ac.lim(:)), -ac.lim(:));
dcmd(4)   = max(min(dcmd(4), ac.thrMax), ac.thrMin);
end

function s = stateOf6(x, ex)
s.r = x(1:3);  s.V = ex.V;  s.C = reshape(x(7:15),3,3);
s.qint = x(19);  s.pint = x(20);  s.nAct = ex.nz;  s.alpha = ex.alpha;  s.w = x(16:18);
C = s.C;
s.phi   = atan2(C(3,2), C(3,3));
s.theta = -asin(max(min(C(3,1),1),-1));
s.psi   = atan2(C(2,1), C(1,1));
vn = C*x(4:6);
s.gamma = asin(max(min(-vn(3)/max(norm(vn),1e-3),1),-1));
s.h = -x(3);
s.upright = C(3,3) > 0;
end

function [al, de, th] = trimLevel(ac, V, h, g0)
% straight and level (theta = alpha): Newton on [alpha, elevator] for the
% body-z force balance and Cm = 0, then thrust from the body-x balance
z = [0.05; 0];
res = @(z) trimRes(ac, z, V, h, g0);
for it = 1:60
    r0 = res(z);
    J = zeros(2);
    for j = 1:2
        dz = zeros(2,1); dz(j) = 1e-5;
        J(:,j) = (res(z + dz) - r0)/1e-5;
    end
    step = -J\r0;
    z = z + max(min(step, 0.05), -0.05);
    if norm(r0) < 1e-10, break; end
end
al = z(1);  de = z(2);
A = ac.aero(al, 0, V, [0;0;0], [de;0;0;0], h, false);
qS = 0.5*A.rho*V^2*A.S;
th = (A.m*g0*sin(al) - A.CX*qS)/(A.StaticThrust*A.lapse);
end

function r = trimRes(ac, z, V, h, g0)
A  = ac.aero(z(1), 0, V, [0;0;0], [z(2);0;0;0], h, false);
qS = 0.5*A.rho*V^2*A.S;
r  = [A.CZ*qS/A.m + g0*cos(z(1)); A.Cm];
end

function C = cnbE(ph, th, ps)
sp = sin(ph); cp = cos(ph); st = sin(th); ct = cos(th); ss = sin(ps); cs = cos(ps);
C = [ct*cs, ct*ss, -st; -cp*ss+sp*st*cs, cp*cs+sp*st*ss, sp*ct; sp*ss+cp*st*cs, -sp*cs+cp*st*ss, cp*ct];
end
