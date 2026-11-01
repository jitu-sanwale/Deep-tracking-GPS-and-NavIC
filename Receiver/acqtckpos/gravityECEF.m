function g = gravityECEF(r)
%Purpose:
%   Acceleration due to gravity (gravitation + centrifugal) in ECEF, WGS-84
%   with the J2 term (Groves, Gravity_ECEF).
%Inputs:
%   r - ECEF position [m]
%Outputs:
%   g - gravity [m/s^2], 3x1
mu = 3.986004418e14; J2 = 1.082627e-3; R0 = 6378137; omgIE = 7.292115e-5;
r  = r(:);
mr = norm(r);
if mr == 0, g = zeros(3,1); return; end
zs  = 5*(r(3)/mr)^2;
gam = -mu/mr^3*(r + 1.5*J2*(R0/mr)^2*[(1-zs)*r(1); (1-zs)*r(2); (3-zs)*r(3)]);
g   = gam + omgIE^2*[r(1); r(2); 0];
end
