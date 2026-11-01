function d = rxTropo(latDeg, heightM, elDeg)
%RXTROPO Receiver tropospheric correction [m].
%
%   d = rxTropo(latDeg, heightM, elDeg)
%
%   Inputs:
%       latDeg  - receiver geodetic latitude [deg]
%       heightM - receiver height above the ellipsoid [m]
%       elDeg   - satellite elevation [deg]
%   Output:
%       d       - slant tropospheric delay [m]
%
%   Uses the same Saastamoinen zenith delay and cosecant mapping as the
%   simulator (saastamoinenTropo.m) with the simulator's standard atmosphere
%   (1013.25 hPa, 288.15 K, 50 % relative humidity; Simulator\initSettings.m),
%   so the tropospheric delay injected into the IF data is removed by the
%   receiver up to the effect of its position estimate. If the simulator
%   settings are changed, change them here as well.
%   A diverged navigation state (non-finite or off-Earth position) gives no
%   correction instead of aborting the run.
%   Written for this package (J. Sanwale et al., 2026), GPL-3.0-or-later.

if ~all(isfinite([latDeg heightM elDeg])) || abs(latDeg) > 90 || abs(heightM) > 1e5
    d = 0;  return;
end
d = saastamoinenTropo(elDeg, latDeg, heightM, 1013.25, 288.15, 50)*299792458;
end
