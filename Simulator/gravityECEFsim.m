function g = gravityECEFsim(r)
%GRAVITYECEFSIM Gravity (gravitation + centrifugal) in ECEF, WGS-84 with J2
%   (Groves, Gravity_ECEF). Same model as the receiver's
%   acqtckpos/gravityECEF.m, kept here so the simulator is self-contained.
mu = 3.986004418e14; J2 = 1.082627e-3; R0 = 6378137; omgIE = 7.292115e-5;
r  = r(:);
mr = norm(r);
zs  = 5*(r(3)/mr)^2;
gam = -mu/mr^3*(r + 1.5*J2*(R0/mr)^2*[(1-zs)*r(1); (1-zs)*r(2); (3-zs)*r(3)]);
g   = gam + omgIE^2*[r(1); r(2); 0];
end
