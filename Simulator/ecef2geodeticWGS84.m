function [latDeg, lonDeg, h] = ecef2geodeticWGS84(xyz)
%ECEF2GEODETICWGS84 Iterative ECEF -> WGS-84 geodetic conversion (vectorised).
%   [latDeg, lonDeg, h] = ecef2geodeticWGS84(xyz), xyz is N x 3 [m].
%   Iterated to 1e-13 rad, i.e. sub-micrometre round trip with geo2cartd(...,5),
%   so the resampled CSV reproduces the spline positions exactly when the
%   simulator converts it back to ECEF.
a  = 6378137;
f  = 1/298.257223563;
e2 = f*(2-f);
x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
p   = hypot(x, y);
lon = atan2(y, x);
lat = atan2(z, p*(1-e2));
for it = 1:20
    N    = a ./ sqrt(1 - e2*sin(lat).^2);
    h    = p./cos(lat) - N;
    latN = atan2(z, p .* (1 - e2*N./(N+h)));
    done = max(abs(latN - lat)) < 1e-13;
    lat  = latN;
    if done, break; end
end
N  = a ./ sqrt(1 - e2*sin(lat).^2);
h  = p./cos(lat) - N;
latDeg = lat*180/pi;
lonDeg = lon*180/pi;
end
