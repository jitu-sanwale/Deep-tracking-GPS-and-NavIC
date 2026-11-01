function [r, v, Cbe] = insNavECEF(tau, rOld, vOld, CbeOld, fib, wib)
%Purpose:
%   One step of the 6-DOF strapdown inertial navigation equations in ECEF
%   (Groves, Principles of GNSS, Inertial and Multisensor Integrated
%   Navigation Systems, 2nd ed., Nav_equations_ECEF: precise attitude update,
%   average-attitude specific-force transformation).
%Inputs:
%   tau     - time interval [s]
%   rOld    - ECEF position [m], 3x1
%   vOld    - ECEF velocity [m/s], 3x1
%   CbeOld  - body-to-ECEF coordinate transformation matrix
%   fib     - average specific force over the interval, body frame [m/s^2]
%   wib     - average angular rate over the interval, body frame [rad/s]
%Outputs:
%   r, v, Cbe - updated position, velocity and attitude

omgIE = 7.292115e-5;
OmgIE = [0 -omgIE 0; omgIE 0 0; 0 0 0];

% attitude update
aIE    = omgIE*tau;
CEarth = [cos(aIE) sin(aIE) 0; -sin(aIE) cos(aIE) 0; 0 0 1];
aIB    = wib(:)*tau;
mag    = norm(aIB);
A      = [0 -aIB(3) aIB(2); aIB(3) 0 -aIB(1); -aIB(2) aIB(1) 0];
if mag > 1e-8
    Cno  = eye(3) + sin(mag)/mag*A + (1-cos(mag))/mag^2*A*A;
    Cave = eye(3) + (1-cos(mag))/mag^2*A + (1-sin(mag)/mag)/mag^2*A*A;
else
    Cno  = eye(3) + A;
    Cave = eye(3);
end
Cbe = CEarth*CbeOld*Cno;

% specific force frame transformation and velocity / position update
CbeAve = CbeOld*Cave - 0.5*OmgIE*CbeOld*tau;
fie = CbeAve*fib(:);
v   = vOld + tau*(fie + gravityECEF(rOld) - 2*OmgIE*vOld);
r   = rOld + 0.5*tau*(v + vOld);
end
