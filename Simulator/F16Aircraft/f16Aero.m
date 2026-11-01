function A = f16Aero(alphar, betar, V, pqr, d, h, withDeriv)
%F16AERO Nonlinear F-16 aerodynamics (NASA TP-1538 tables, Stevens & Lewis).
%   A = f16Aero(alphar, betar, V, pqr, d, h, withDeriv)
%
%   MATLAB port of the aerodynamic build-up of F16Model.jl
%   (rfStengel_AircraftModel\F16Model.jl-master, MIT License, Intelligent
%   Systems Research Laboratory 2020; see LICENSE_F16Model_jl.txt), which
%   follows the UMN/Stevens & Lewis "nlplant" F-16 model. Data:
%   F16AeroData.h5 (copied from that package), multilinear interpolation on
%   the table grids as in the package, clamped at the grid edges
%   (alpha -20..90 deg / LEF tables -20..45 deg, beta -30..30 deg).
%   Differences from F16Model.jl:
%     - Z-force pitch damping uses delta_Czq_lef (the Julia code uses
%       delta_Cz_lef there; the original nlplant.c uses delta_Czq_lef).
%     - The leading-edge flap follows the standard F-16 automatic schedule
%       lef = 1.38*alpha - 9.05*qbar/ps + 1.45 deg, limited to 0..25 deg.
%     - SI units in and out.
%
%   alphar, betar [rad]; V airspeed [m/s]; pqr body rates [rad/s];
%   d = [elevator aileron rudder thrust] ([rad rad rad], thrust as a
%   fraction of 19000 lbf); h altitude [m]
%   withDeriv - also return the control effectiveness (finite difference for
%   the elevator, exact for aileron/rudder, which enter linearly)
%   A - CX CY CZ Cl Cm Cn (body axes), Thrust [N], rho, mass properties
%       (m Ixx Iyy Izz Ixz S b cBar), StaticThrust, lapse, and with
%       withDeriv: dCm_de dCl_da dCl_dr dCn_da dCn_dr (per rad)

persistent T
if isempty(T), T = loadTables(); end
if nargin < 7, withDeriv = false; end

r2d = 180/pi;
al = min(max(alphar*r2d, -20), 90);  be = min(max(betar*r2d, -30), 30);
el = d(1)*r2d;  ail = d(2)*r2d;  rud = d(3)*r2d;

% atmosphere (F16Model.jl atmos, English units) -> qbar/ps for the LEF schedule
hft = h/0.3048;  tfac = 1 - 0.703e-5*hft;  temp = 519*tfac;  if hft >= 35000, temp = 390; end
rhoE = 2.377e-3*tfac^4.14;  ps = 1715*rhoE*temp;
vft = max(V, 0.1)/0.3048;   qbarE = 0.5*rhoE*vft^2;
lef = min(max(1.38*al - 9.05*qbarE/ps + 1.45, 0), 25);

c = coeffs(T, al, be, el, ail, rud, lef, vft, pqr);
A.CX = c(1); A.CY = c(2); A.CZ = c(3); A.Cl = c(4); A.Cm = c(5); A.Cn = c(6);
A.lef = lef;
if withDeriv
    de = 1;                                              % deg
    cp = coeffs(T, al, be, min(el+de,25), ail, rud, lef, vft, pqr);
    cm = coeffs(T, al, be, max(el-de,-25), ail, rud, lef, vft, pqr);
    A.dCm_de = (cp(5) - cm(5))/((min(el+de,25) - max(el-de,-25))/r2d);
    A.dCl_da = c(7)*r2d;  A.dCn_da = c(8)*r2d;           % per rad
    A.dCl_dr = c(9)*r2d;  A.dCn_dr = c(10)*r2d;
end

A.rho = rhoE*515.379;                                    % slug/ft^3 -> kg/m^3
A.StaticThrust = 19000*4.44822;  A.lapse = 1;            % thrust input as in the package
A.Thrust = d(4)*A.StaticThrust;
A.m   = 636.94*14.5939;                                  % kg
A.S   = 300*0.3048^2;  A.b = 30*0.3048;  A.cBar = 11.32*0.3048;
kI = 14.5939*0.3048^2;                                   % slug ft^2 -> kg m^2
A.Ixx = 9496*kI;  A.Iyy = 55814*kI;  A.Izz = 63100*kI;  A.Ixz = 982*kI;
end

%% ------------------------------------------------------------------------
function c = coeffs(T, al, be, el, ail, rud, lef, vt, pqr)
% Returns [CX CY CZ Cl Cm Cn dCl/dail dCn/dail dCl/drud dCn/drud] (per deg)
P = pqr(1); Q = pqr(2); R = pqr(3);
B = 30; cbar = 11.32; xcgr = 0.35; xcg = 0.30;
al2 = min(max(al, -20), 45);                             % LEF tables stop at 45 deg
dail = ail/21.5;  drud = rud/30;  dlef = 1 - lef/25;

% grouped table look-ups (multi-valued griddedInterpolants, see loadTables)
v  = T.G3m([al be el; al be 0]);                         % Cx Cz Cm at el and at 0
w  = T.G3l([al be el; al be 0]);                         % Cl Cn at el and at 0
u2 = T.G2a(al, be);                                      % Cy Cy_r30 Cn_r30 Cl_r30 Cy_a20 Cn_a20 Cl_a20
l2 = T.G2l(al2, be);                                     % Cx Cy Cz Cl Cm Cn _lef, Cy Cn Cl _a20_lef
d1 = T.G1a(al);                                          % Cxq Cyp Czq Cmq Cyr Cnr Cnp Clp Clr dCnbeta dClbeta dCm
e1 = T.G1l(al2);                                         % delta Cxq Cyr Cyp Czq Clr Clp Cmq Cnr Cnp _lef
eta = T.eta_el(el);
Cx = v(1,1); Cz = v(1,2); Cm = v(1,3); Cx0 = v(2,1); Cz0 = v(2,2); Cm0 = v(2,3);
Cl = w(1,1); Cn = w(1,2); Cl0 = w(2,1); Cn0 = w(2,2);
Cy = u2(1);

dCx_lef = l2(1) - Cx0;  dCy_lef = l2(2) - Cy;  dCz_lef = l2(3) - Cz0;
dCl_lef = l2(4) - Cl0;  dCm_lef = l2(5) - Cm0;  dCn_lef = l2(6) - Cn0;

Cxq = d1(1); Cyp = d1(2); Czq = d1(3); Cmq = d1(4); Cyr = d1(5); Cnr = d1(6); Cnp = d1(7); Clp = d1(8); Clr = d1(9);
dCxq = e1(1); dCyr = e1(2); dCyp = e1(3); dCzq = e1(4); dClr = e1(5); dClp = e1(6); dCmq = e1(7); dCnr = e1(8); dCnp = e1(9);

dCy_r30 = u2(2) - Cy;   dCn_r30 = u2(3) - Cn0;  dCl_r30 = u2(4) - Cl0;
dCy_a20 = u2(5) - Cy;   dCy_a20_lef = l2(7) - l2(2) - dCy_a20;
dCn_a20 = u2(6) - Cn0;  dCn_a20_lef = l2(8) - l2(6) - dCn_a20;
dCl_a20 = u2(7) - Cl0;  dCl_a20_lef = l2(9) - l2(4) - dCl_a20;
dCnbeta = d1(10);  dClbeta = d1(11);  dCm = d1(12);

Cx_tot = Cx + dCx_lef*dlef + (cbar/(2*vt))*(Cxq + dCxq*dlef)*Q;
Cz_tot = Cz + dCz_lef*dlef + (cbar/(2*vt))*(Czq + dCzq*dlef)*Q;
Cm_tot = Cm*eta + Cz_tot*(xcgr - xcg) + dCm_lef*dlef + (cbar/(2*vt))*(Cmq + dCmq*dlef)*Q + dCm;
dYdail = dCy_a20 + dCy_a20_lef*dlef;
Cy_tot = Cy + dCy_lef*dlef + dYdail*dail + dCy_r30*drud + (B/(2*vt))*(Cyr + dCyr*dlef)*R + (B/(2*vt))*(Cyp + dCyp*dlef)*P;
dNdail = dCn_a20 + dCn_a20_lef*dlef;
Cn_tot = Cn + dCn_lef*dlef - Cy_tot*(xcgr - xcg)*(cbar/B) + dNdail*dail + dCn_r30*drud ...
         + (B/(2*vt))*(Cnr + dCnr*dlef)*R + (B/(2*vt))*(Cnp + dCnp*dlef)*P + dCnbeta*be;
dLdail = dCl_a20 + dCl_a20_lef*dlef;
Cl_tot = Cl + dCl_lef*dlef + dLdail*dail + dCl_r30*drud + (B/(2*vt))*(Clr + dClr*dlef)*R ...
         + (B/(2*vt))*(Clp + dClp*dlef)*P + dClbeta*be;
% control effectiveness per degree (aileron normalised by 21.5, rudder by 30);
% the -Cy*(xcgr-xcg)*cbar/B coupling of Cn is included
dCn_dail = dNdail/21.5 - (dYdail/21.5)*(xcgr - xcg)*(cbar/B);
dCn_drud = dCn_r30/30  - (dCy_r30/30)*(xcgr - xcg)*(cbar/B);
c = [Cx_tot Cy_tot Cz_tot Cl_tot Cm_tot Cn_tot dLdail/21.5 dCn_dail dCl_r30/30 dCn_drud];
end

function T = loadTables()
f = fullfile(fileparts(mfilename('fullpath')), 'F16AeroData.h5');
rd = @(n) double(h5read(f, ['/' n]));
a1 = rd('alpha1'); a2 = rd('alpha2'); b1 = rd('beta1'); d1 = rd('dh1'); d2 = rd('dh2');
% tables sharing a grid are stacked along a trailing dimension so one
% griddedInterpolant call returns all of them (same multilinear values)
s3 = @(names, a, b, e) stackTables(rd, names, [numel(a) numel(b) numel(e)]);
s2 = @(names, a, b) stackTables(rd, names, [numel(a) numel(b)]);
s1 = @(names) stackTables(rd, names, []);
T.G3m = griddedInterpolant({a1, b1, d1}, s3({'_Cx','_Cz','_Cm'}, a1, b1, d1), 'linear', 'nearest');
T.G3l = griddedInterpolant({a1, b1, d2}, s3({'_Cl','_Cn'}, a1, b1, d2), 'linear', 'nearest');
T.G2a = griddedInterpolant({a1, b1}, s2({'_Cy','_Cy_r30','_Cn_r30','_Cl_r30','_Cy_a20','_Cn_a20','_Cl_a20'}, a1, b1), 'linear', 'nearest');
T.G2l = griddedInterpolant({a2, b1}, s2({'_Cx_lef','_Cy_lef','_Cz_lef','_Cl_lef','_Cm_lef','_Cn_lef', ...
                           '_Cy_a20_lef','_Cn_a20_lef','_Cl_a20_lef'}, a2, b1), 'linear', 'nearest');
T.G1a = griddedInterpolant(a1, s1({'_Cxq','_Cyp','_Czq','_Cmq','_Cyr','_Cnr','_Cnp','_Clp','_Clr', ...
                           '_deltaCnbeta','_deltaClbeta','_deltaCm'}), 'linear', 'nearest');
T.G1l = griddedInterpolant(a2, s1({'_deltaCxq_lef','_deltaCyr_lef','_deltaCyp_lef','_deltaCzq_lef','_deltaClr_lef', ...
                           '_deltaClp_lef','_deltaCmq_lef','_deltaCnr_lef','_deltaCnp_lef'}), 'linear', 'nearest');
T.eta_el = griddedInterpolant(d1, rd('_eta_el'), 'linear', 'nearest');
T.alpha1 = a1;  T.alpha2 = a2;
end

function V = stackTables(rd, names, sz)
% tables reshaped to the grid size sz and stacked along the next dimension
% (sz = [] for 1-D tables: columns of an n x k matrix)
parts = cell(1, numel(names));
for k = 1:numel(names)
    x = rd(names{k});
    if isempty(sz), parts{k} = x(:); else, parts{k} = reshape(x, sz); end
end
V = cat(max(numel(sz), 1) + 1, parts{:});
if isempty(sz), V = cat(2, parts{:}); end
end
