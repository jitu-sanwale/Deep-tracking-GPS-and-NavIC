function tropoDelaySec = saastamoinenTropo(elDeg, latDeg, heightM, P0, T0, RH)
%SAASTAMOINENTROPO Slant tropospheric delay (seconds) via the
%   Saastamoinen zenith-total-delay model.
%
%   tropoDelaySec = saastamoinenTropo(elDeg, latDeg, heightM, P0, T0, RH)
%
%   Reference: J. Saastamoinen, "Contributions to the theory of
%   atmospheric refraction: Part II," Bull. Geod., vol. 107, 1973,
%   doi:10.1007/BF02522083.
%
%   Inputs:
%       elDeg    - satellite elevation [deg]
%       latDeg   - receiver geodetic latitude [deg]
%       heightM  - receiver height above the ellipsoid [m]
%       P0       - surface pressure [mbar/hPa]
%       T0       - surface temperature [K]
%       RH       - relative humidity [%]
%
%   Output:
%       tropoDelaySec - slant tropospheric group delay [s]
%
%   Note: this is intentionally a DIFFERENT model from the receiver's own
%   tropospheric correction (Common/tropo.m, a Goad & Goodman model), the
%   same way a real receiver's tropo model never exactly matches the true
%   atmosphere - some residual error after correction is expected and
%   physically realistic.

c = 299792458;

% Saturation water vapor pressure (Magnus formula, T in Celsius).
Tc = T0 - 273.15;
e0 = (RH/100) * 6.11 * 10^(7.5*Tc/(237.3+Tc));

latRad = deg2rad(latDeg);
hKm = heightM/1000;

% Zenith total delay [m].
ztdMeters = 0.002277/(1 - 0.00266*cos(2*latRad) - 0.00028*hKm) * (P0 + (1255/T0 + 0.05)*e0);

% Simple cosecant obliquity mapping (elevation floored to avoid the
% singularity at the horizon). More advanced mapping functions such as
% Niell/VMF3/GPT3 require additional meteorological grids not available
% here.
elClamped = max(elDeg, 5);
slantMeters = ztdMeters / sind(elClamped);

tropoDelaySec = slantMeters / c;
