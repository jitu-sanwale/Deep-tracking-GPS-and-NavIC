function llh = xyz2llh(xyz)
%XYZ2LLH WGS-84 geodetic latitude, longitude and height from ECEF position.
%
%   llh = xyz2llh(xyz)
%
%   Input:
%       xyz - ECEF position [m] (row or column)
%   Output:
%       llh - 1x3 row: [latitude (rad), longitude (rad), ellipsoidal height (m)]
%
%   Fixed-point iteration on the geodetic latitude (e.g. Groves, Principles
%   of GNSS, Inertial, and Multisensor Integrated Navigation Systems, 2nd ed.,
%   2013); converges to machine precision in a few steps for points
%   near the Earth. Written for this package (J. Sanwale et al., 2026),
%   GPL-3.0-or-later.

a  = 6378137.0;                 % WGS-84 semi-major axis [m]
f  = 1/298.257223563;           % WGS-84 flattening
e2 = f*(2 - f);                 % first eccentricity squared

x = xyz(1);  y = xyz(2);  z = xyz(3);
p   = hypot(x, y);
lon = atan2(y, x);
lat = atan2(z, p*(1 - e2));
for k = 1:20
    N   = a/sqrt(1 - e2*sin(lat)^2);
    h   = p/cos(lat) - N;
    latNew = atan2(z, p*(1 - e2*N/(N + h)));
    converged = abs(latNew - lat) < 1e-15;
    lat = latNew;
    if converged, break; end
end
N = a/sqrt(1 - e2*sin(lat)^2);
if abs(lat) < pi/4
    h = p/cos(lat) - N;
else
    h = z/sin(lat) - N*(1 - e2);   % better conditioned near the poles
end
llh = [lat, lon, h];
end
