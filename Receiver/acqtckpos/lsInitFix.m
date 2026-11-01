function [pos, clkBias, ok, isb] = lsInitFix(pr, svxyz, isNavIC)
%Purpose:
%   Iterative least-squares single-point fix (position + clock bias, both in
%   meters) used to initialise the EKF. Earth-rotation correction included.
%   With NavIC channels (isNavIC true) and at least five satellites, the
%   GPS-NavIC inter-system bias isb (m) is a fifth unknown:
%       pr_i = |r_i - x| + b + isNavIC_i * isb
%Inputs:
%   pr      - pseudoranges corrected for SV clock/TGD (and, for NavIC, the
%             IRNSS-GPS time offset of message type 9) (m), 1xN
%   svxyz   - SV positions at transmit time, Nx3 (ECEF, m)
%   isNavIC - optional 1xN flags (default: all GPS)
n  = numel(pr);
if nargin < 3 || isempty(isNavIC), isNavIC = zeros(1,n); end
isNavIC = double(isNavIC(:)');
useIsb = any(isNavIC) && n >= 5;
nx = 4 + useIsb;
x  = zeros(nx,1);
ok = false;
for it = 1:30
    H  = zeros(n,nx);
    dz = zeros(n,1);
    for k = 1:n
        svr     = erotcorr(svxyz(k,:), pr(k) - x(4));
        rho     = norm(svr(:) - x(1:3));
        H(k,1:4) = [-(svr(:)' - x(1:3)')/rho, 1];
        bk      = x(4);
        if useIsb, H(k,5) = isNavIC(k); bk = bk + isNavIC(k)*x(5); end
        dz(k)   = pr(k) - rho - bk;
    end
    dx = H \ dz;
    x  = x + dx;
    if norm(dx) < 1e-4, ok = true; break; end
end
pos     = x(1:3)';
clkBias = x(4);
isb     = 0;
if useIsb, isb = x(5); end
end
