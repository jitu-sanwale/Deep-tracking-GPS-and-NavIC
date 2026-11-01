function [azBody, elBody] = applyAttitudeToBodyFrame(azDeg, elDeg, rollDeg, pitchDeg, headingDeg)
%APPLYATTITUDETOBODYFRAME Rotate a satellite direction from the local
%   level frame into the receive antenna's body frame, given the
%   platform's attitude.
%
%   [azBody, elBody] = applyAttitudeToBodyFrame(azDeg, elDeg, rollDeg, pitchDeg, headingDeg)
%
%   Inputs azDeg/elDeg are local-level (topocent.m's convention: azimuth
%   clockwise from North, elevation above the local horizon). Attitude is
%   the standard aerospace 3-2-1 Euler set relative to local North-East-
%   Down: heading (yaw) clockwise from North, pitch positive nose-up, roll
%   positive right-wing-down. Body axes are x forward, y right, z down; the
%   antenna boresight is body -z (the top of the airframe). Outputs are the
%   direction's azimuth (clockwise from the nose, about body -z) and
%   elevation above the antenna's horizon - at zero attitude they equal
%   the inputs.
%
%   Replaces a version that applied the rotations to a [North East Up]
%   vector as if it were right-handed: that flipped the sign of pitch and
%   roll (nose-up 30 deg put a satellite straight ahead on the horizon at
%   +30 deg body elevation instead of -30 deg). Heading-only and zero
%   attitude were unaffected.

vN = [cosd(elDeg)*cosd(azDeg); cosd(elDeg)*sind(azDeg); -sind(elDeg)];   % NED unit vector

sp = sind(rollDeg);    cp = cosd(rollDeg);
st = sind(pitchDeg);   ct = cosd(pitchDeg);
ss = sind(headingDeg); cs = cosd(headingDeg);
Cnb = [ct*cs,               ct*ss,               -st;
       -cp*ss + sp*st*cs,   cp*cs + sp*st*ss,    sp*ct;
       sp*ss + cp*st*cs,   -sp*cs + cp*st*ss,    cp*ct];                 % NED -> body

vB = Cnb*vN;
elBody = asind(max(min(-vB(3), 1), -1));
azBody = mod(atan2d(vB(2), vB(1)), 360);
end
