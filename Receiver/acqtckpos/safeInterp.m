function yi = safeInterp(x, y, xi)
%Purpose:
%   Linear interpolation between two (sample-index, value) points that falls
%   back to the newest value when the points are degenerate or xi lies far
%   outside them (happens when a channel loses lock and its code NCO runs away;
%   plain interp1 extrapolation then produced values of -1e7 samples).
dx = x(2) - x(1);
if ~all(isfinite([x(:); y(:)])) || dx < 1
    yi = y(end); return;
end
f = (xi - x(1)) / dx;
if f < -0.5 || f > 1.5
    yi = y(end);
else
    yi = y(1) + f*(y(2) - y(1));
end
end
