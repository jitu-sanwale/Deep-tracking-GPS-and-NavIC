function enu = xyz2enu(xyz, orgxyz)
%XYZ2ENU Local east-north-up coordinates of an ECEF point about an origin.
%
%   enu = xyz2enu(xyz, orgxyz)
%
%   Inputs:
%       xyz    - ECEF position of the point [m] (row or column)
%       orgxyz - ECEF position of the local origin [m] (row or column)
%   Output:
%       enu    - 3x1 column: [east; north; up] relative to the origin [m],
%                on the WGS-84 geodetic tangent plane at the origin
%
%   Written for this package (J. Sanwale et al., 2026), GPL-3.0-or-later.

if nargin < 2, error('xyz2enu: point and origin are both required'); end
d = xyz(:) - orgxyz(:);
llh = xyz2llh(orgxyz);
sl = sin(llh(1));  cl = cos(llh(1));
so = sin(llh(2));  co = cos(llh(2));
R = [ -so,     co,    0;
      -sl*co, -sl*so, cl;
       cl*co,  cl*so, sl ];
enu = R*d;
end
