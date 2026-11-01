function ionodel = ionocorr(gpsTow, svXyz, usrXyz)
%IONOCORR Broadcast (Klobuchar) single-frequency ionospheric delay, L1.
%
%   ionodel = ionocorr(gpsTow, svXyz, usrXyz)
%
%   Inputs:
%       gpsTow - GPS time of week at which the delay is evaluated [s]
%       svXyz  - satellite position, ECEF [m] (row or column)
%       usrXyz - user position, ECEF [m] (row or column)
%   Global:
%       ALPHA, BETA - 1x4 broadcast ionospheric coefficients
%   Output:
%       ionodel - L1 ionospheric group delay [m]; the caller scales it to
%                 other frequencies (NavIC L5: gamma = (f_L1/f_L5)^2)
%
%   Model: J. A. Klobuchar, IEEE Trans. Aerosp. Electron. Syst., AES-23(3),
%   1987; IS-GPS-200, 20.3.3.5.2.5. Angles inside the model are in
%   semicircles. The geomagnetic-latitude term uses 1.617 (IS-GPS-200), as in
%   the simulator (klobuchar.m), so the receiver removes the ionospheric delay
%   the simulator injects, given the same broadcast coefficients.
%   Written for this package (J. Sanwale et al., 2026), GPL-3.0-or-later.

global ALPHA BETA

usr = usrXyz(:);
los = svXyz(:) - usr;

% user geodetic latitude and longitude (WGS-84)
[lat, lon] = ecefToGeodetic(usr);

% satellite elevation and azimuth in the local east-north-up frame
sl = sin(lat); cl = cos(lat); so = sin(lon); co = cos(lon);
e = -so*los(1) + co*los(2);
n = -sl*co*los(1) - sl*so*los(2) + cl*los(3);
u =  cl*co*los(1) + cl*so*los(2) + sl*los(3);
el = atan2(u, hypot(e, n));
az = atan2(e, n);

E = el/pi;                                   % elevation [semicircles]
psi = 0.0137/(E + 0.11) - 0.022;             % earth-centred angle [semicircles]
phiI = min(max(lat/pi + psi*cos(az), -0.416), 0.416);   % pierce-point latitude
lamI = lon/pi + psi*sin(az)/cos(phiI*pi);    % pierce-point longitude
phiM = phiI + 0.064*cos((lamI - 1.617)*pi);  % geomagnetic latitude

t = mod(4.32e4*lamI + gpsTow, 86400);        % local time at the pierce point [s]
F = 1 + 16*(0.53 - E)^3;                     % obliquity factor

AMP = max(ALPHA(1) + ALPHA(2)*phiM + ALPHA(3)*phiM^2 + ALPHA(4)*phiM^3, 0);
PER = max(BETA(1) + BETA(2)*phiM + BETA(3)*phiM^2 + BETA(4)*phiM^3, 72000);
x = 2*pi*(t - 50400)/PER;

if abs(x) < 1.57
    tIono = F*(5e-9 + AMP*(1 - x^2/2 + x^4/24));
else
    tIono = F*5e-9;
end
ionodel = tIono*299792458;
end

function [lat, lon] = ecefToGeodetic(r)
% WGS-84 geodetic latitude and longitude [rad] by fixed-point iteration on
% the latitude (converges to machine precision in a few steps near Earth).
a = 6378137.0;  f = 1/298.257223563;  e2 = f*(2 - f);
p = hypot(r(1), r(2));
lon = atan2(r(2), r(1));
lat = atan2(r(3), p*(1 - e2));
for k = 1:10
    N = a/sqrt(1 - e2*sin(lat)^2);
    h = p/cos(lat) - N;
    latNew = atan2(r(3), p*(1 - e2*N/(N + h)));
    if abs(latNew - lat) < 1e-14, lat = latNew; break; end
    lat = latNew;
end
end
