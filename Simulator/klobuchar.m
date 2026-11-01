function ionoDelaySec = klobuchar(alpha, beta, gpsTimeSec, userLatDeg, userLonDeg, elDeg, azDeg)
%KLOBUCHAR Ionospheric group delay (L1, seconds) via the broadcast
%   Klobuchar model.
%
%   ionoDelaySec = klobuchar(alpha, beta, gpsTimeSec, userLatDeg, userLonDeg, elDeg, azDeg)
%
%   Reference: J.A. Klobuchar, "Ionospheric Time-Delay Algorithm for
%   Single-Frequency GPS Users," IEEE Trans. Aerosp. Electron. Syst.,
%   AES-23(3), 1987, doi:10.1109/TAES.1987.310829 (also IS-GPS-200,
%   20.3.3.5.2.5). Identical implementation is used on both the
%   simulator (to inject the delay) and the receiver (to correct it),
%   given the same alpha/beta coefficients.
%
%   Inputs:
%       alpha, beta        - 1x4 broadcast coefficient vectors (from the
%                            nav message / RINEX header ION ALPHA / ION
%                            BETA lines)
%       gpsTimeSec          - GPS time of week [s]
%       userLatDeg/LonDeg  - user geodetic position [deg]
%       elDeg, azDeg        - satellite elevation/azimuth as seen from
%                            the user [deg]
%
%   Output:
%       ionoDelaySec - single-frequency (L1) ionospheric group delay [s]

% The model's native unit is semicircles (degrees/180).
lat = userLatDeg/180;
lon = userLonDeg/180;
el  = elDeg/180;
az  = deg2rad(azDeg);

% Earth-centered angle between the user and the ionospheric pierce point.
psi = 0.0137/(el+0.11) - 0.022;

% Subionospheric point geomagnetic latitude/longitude.
phi_i = lat + psi*cos(az);
phi_i = max(min(phi_i, 0.416), -0.416);
lambda_i = lon + psi*sin(az)/cos(phi_i*pi);
phi_m = phi_i + 0.064*cos((lambda_i - 1.617)*pi);

% Local time at the subionospheric point.
t = mod(4.32e4*lambda_i + gpsTimeSec, 86400);

% Obliquity (slant) factor.
F = 1.0 + 16.0*(0.53 - el)^3;

% Amplitude and period of the (cosine-modeled) diurnal delay term.
AMP = alpha(1) + alpha(2)*phi_m + alpha(3)*phi_m^2 + alpha(4)*phi_m^3;
AMP = max(AMP, 0);
PER = beta(1) + beta(2)*phi_m + beta(3)*phi_m^2 + beta(4)*phi_m^3;
PER = max(PER, 72000);

x = 2*pi*(t - 50400)/PER;

if abs(x) < 1.57
    ionoDelaySec = F*(5e-9 + AMP*(1 - x^2/2 + x^4/24));
else
    ionoDelaySec = F*5e-9;
end
