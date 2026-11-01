function svxyzr = erotcorr(svxyz, pr)
%EROTCORR Earth-rotation (Sagnac) correction of a satellite position.
%
%   svxyzr = erotcorr(svxyz, pr)
%
%   Inputs:
%       svxyz - satellite ECEF position at transmission, in the ECEF frame
%               of the transmission epoch [m] (row or column)
%       pr    - signal travel distance [m], i.e. the range or the pseudorange
%               corrected for the satellite clock; the travel time is pr/c
%   Output:
%       svxyzr - 3x1 column: the same position expressed in the ECEF frame of
%                the reception epoch, rotated about z by omega_ie*pr/c
%
%   Rotation over the signal travel time as in IS-GPS-200 and Groves (2013).
%   Written for this package (J. Sanwale et al., 2026), GPL-3.0-or-later.

omegaE = 7.2921151467e-5;       % WGS-84 Earth rotation rate [rad/s]
c      = 299792458;             % speed of light [m/s]
theta  = omegaE*pr/c;
R3 = [ cos(theta), sin(theta), 0;
      -sin(theta), cos(theta), 0;
       0,          0,          1 ];
svxyzr = R3*svxyz(:);
end
