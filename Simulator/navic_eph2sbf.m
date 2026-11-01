function payload = navic_eph2sbf(eph, ionoutc, igt)
%NAVIC_EPH2SBF Build the static part of the four 292-bit NavIC L5 SPS
%subframes (before TOWC/CRC/tail/FEC/interleaving) from one satellite's
%ephemeris, per IRNSS SPS ICD v1.1 Table 11 (subframe 1), Table 12
%(subframe 2) and section 6.1.3 (subframe 3 & 4, message type 11).
%
% Bit numbering below matches NavIC_Software_Receiver/include/ephemeris.m
% exactly (1-indexed, MSB-first), which is what this must decode into.
%
%   payload  - 4x292 logical/double matrix, rows = [SF1;SF2;SF3;SF4].
%              Bits 9:25 (TOWC) and 263:292 (CRC+tail) are left as 0 -
%              navic_generateNavMsg.m fills those in per-instant.

gpsPi = 3.1415926535898;

payload = zeros(4, 292);

subframeIdBits = [0 0; 0 1; 1 0; 1 1]; % SF1..SF4 -> Table 10
for k = 1:4
    payload(k, 28:29) = subframeIdBits(k, :);
end

%% ---------------- Subframe 1 : Table 11 (232 data bits, 31:262) --------
ura = min(max(round(eph.svaccur), 0), 15);
l5flag = double(eph.svhealth ~= 0);
sflag  = double(eph.svhealth ~= 0);

sf1 = [ ...
    ubits_(eph.weekNrm, 10), ...
    sbits_(eph.af0 / 2^-31, 22), ...
    sbits_(eph.af1 / 2^-43, 16), ...
    sbits_(eph.af2 / 2^-55, 8), ...
    ubits_(ura, 4), ...
    ubits_(eph.toc / 16, 16), ...
    sbits_(eph.tgd / 2^-31, 8), ...
    sbits_(eph.deltan / (2^-41 * gpsPi), 22), ...
    ubits_(eph.IODEC, 8), ...
    zeros(1, 10), ...            % reserved
    l5flag, sflag, ...
    sbits_(eph.cuc / 2^-28, 15), ...
    sbits_(eph.cus / 2^-28, 15), ...
    sbits_(eph.cic / 2^-28, 15), ...
    sbits_(eph.cis / 2^-28, 15), ...
    sbits_(eph.crc / 2^-4, 15), ...
    sbits_(eph.crs / 2^-4, 15), ...
    sbits_(eph.idot / (2^-43 * gpsPi), 14), ...
    zeros(1, 2) ...               % spare
    ];
assert(numel(sf1) == 232, 'navic_eph2sbf: SF1 data field width mismatch');
payload(1, 31:262) = sf1;

%% ---------------- Subframe 2 : Table 12 (232 data bits, 31:262) --------
sf2 = [ ...
    sbits_(eph.M0 / (2^-31 * gpsPi), 32), ...
    ubits_(eph.toe / 16, 16), ...
    ubits_(eph.ecc / 2^-33, 32), ...
    ubits_(eph.sqrta / 2^-19, 32), ...
    sbits_(eph.Omega0 / (2^-31 * gpsPi), 32), ...
    sbits_(eph.omega / (2^-31 * gpsPi), 32), ...
    sbits_(eph.Omegadot / (2^-41 * gpsPi), 22), ...
    sbits_(eph.i0 / (2^-31 * gpsPi), 32), ...
    zeros(1, 2) ...                % spare
    ];
assert(numel(sf2) == 232, 'navic_eph2sbf: SF2 data field width mismatch');
payload(2, 31:262) = sf2;

%% ---------------- Subframe 3 & 4 : message type 11 (EOP + iono) --------
% Message data occupies absolute bits 37:256 (220 bits). Only the
% ionospheric alpha/beta sub-fields (absolute bits 175:238) are ones the
% receiver's ephemeris.m actually decodes; the EOP portion is left as 0.
data220 = zeros(1, 220);
ionoFields = { ...
    ionoutc.alpha0, 2^-30,        175; ...
    ionoutc.alpha1, 2^-27/gpsPi,  183; ...
    ionoutc.alpha2, 2^-24/gpsPi^2,191; ...
    ionoutc.alpha3, 2^-24/gpsPi^3,199; ...
    ionoutc.beta0,  2^11,         207; ...
    ionoutc.beta1,  2^14/gpsPi,   215; ...
    ionoutc.beta2,  2^16/gpsPi^2, 223; ...
    ionoutc.beta3,  2^16/gpsPi^3, 231  ...
    };
for r = 1:size(ionoFields, 1)
    val = ionoFields{r, 1};
    scale = ionoFields{r, 2};
    absStart = ionoFields{r, 3};
    relStart = absStart - 36; % offset into data220 (data220(1) == abs bit 37)
    data220(relStart : relStart + 7) = sbits_(val / scale, 8);
end

for k = 3:4
    payload(k, 31:36) = ubits_(11, 6);        % Message ID = 11
    payload(k, 37:256) = data220;
    payload(k, 257:262) = ubits_(eph.PRN, 6); % PRN ID
end

%% ---------------- Subframe 3 : message type 9 (IRNSS - GPS time offset) ---
% Optional third argument igt (struct A0 [s], A1 [s/s], A2 [s/s^2], tot [s of
% week], WNot [IRNSS week]): subframe 3 then carries message type 9 (ICD
% Table 16 / Figure 17) - IRNSS-UTC fields zero, IRNSS-GPS fields A0 (bits
% 129-144, 2^-35 s), A1 (145-157, 2^-51), A2 (158-164, 2^-68), tot (165-180,
% x16 s), WNot (181-190), GNSS ID = 0 for GPS (191-193), spare (194-256).
% The user recovers GPS time as t_GPS = t_IRN - [A0 + A1 dt + A2 dt^2],
% dt = t_IRN - tot + 604800 (WN - WNot) (ICD Appendix F.2).
if nargin >= 3 && ~isempty(igt)
    mt9 = zeros(1, 292);
    mt9(129:144) = sbits_(igt.A0 / 2^-35, 16);
    mt9(145:157) = sbits_(igt.A1 / 2^-51, 13);
    mt9(158:164) = sbits_(igt.A2 / 2^-68, 7);
    mt9(165:180) = ubits_(igt.tot / 16, 16);
    mt9(181:190) = ubits_(mod(igt.WNot, 1024), 10);
    mt9(191:193) = 0;                          % GNSS ID: GPS
    payload(3, 31:36)   = ubits_(9, 6);        % Message ID = 9
    payload(3, 37:256)  = mt9(37:256);
end

end

%% ------------------------------------------------------------------ %%
function bits = ubits_(u, n)
% Unsigned integer -> n-bit row vector, MSB first. Saturates (rather than
% wrapping) if the value doesn't fit, so an out-of-range real-world
% coefficient degrades gracefully instead of aliasing to a wild value.
u = min(max(round(u), 0), 2^n - 1);
bits = double(bitget(uint64(u), n:-1:1));
end

function bits = sbits_(v, n)
% Two's-complement signed integer -> n-bit row vector, MSB first. Saturates
% to the representable range [-2^(n-1), 2^(n-1)-1] before encoding.
v = min(max(round(v), -2^(n - 1)), 2^(n - 1) - 1);
bits = ubits_(mod(v, 2^n), n);
end
