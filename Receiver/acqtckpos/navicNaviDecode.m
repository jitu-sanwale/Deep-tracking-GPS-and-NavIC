function [e, nav1, sfb1] = navicNaviDecode(Pi, prn)
%NAVICNAVIDECODE Decode the NavIC L5 SPS navigation message of one channel
%   from its 1-ms in-phase prompt correlations (IRNSS SPS ICD v1.1, Ch. 5-6).
%   [e, nav1, sfb1] = navicNaviDecode(Pi, prn)
%     Pi   - P_i of the channel, one value per ms (conventional tracking)
%     e    - ephemeris struct in the receiver's GPS field convention
%            (naviDecode.m / svPosVel.m), plus the IRNSS-GPS time offset of
%            message type 9 (igtA0 A1 A2 tot WNot) and isNavIC = true
%     nav1 - ms (0-based) of the first symbol boundary
%     sfb1 - symbols (0-based, from nav1) to the start of the first subframe
%            that passed its CRC; e.TOW1 is the IRNSS time of that start
%   Steps: symbol synchronisation (sign-change histogram mod 20 ms), 50-sps
%   symbols, EB90 sync search (either polarity, confirmed 600 symbols
%   later), 73x8 block de-interleaving, Viterbi decoding of the rate-1/2 K=7
%   (171,133) code, CRC-24Q, then Table 11/12 (subframes 1, 2) and message
%   type 9 (subframe 3/4). TOW of the subframe start = 12*TOWC - 12.
%   TGD is returned scaled by (fS/fL5)^2, so that svPosVel's "clock - TGD"
%   is the L5 single-frequency correction of IrnssClockCorrection.m.

e = struct('isNavIC', true, 'updateflag', 0);
nav1 = NaN; sfb1 = NaN;
Pi = Pi(:)';
s  = sign(Pi); s(s == 0) = 1;
tr = find(s(2:end) ~= s(1:end-1)) + 1;          % ms on which a new symbol starts
tr = tr(tr > 500);
if numel(tr) < 5, fprintf('NavIC PRN %2d: no symbol transitions\n', prn); return; end
h  = accumarray(mod(tr(:)-1,20)+1, 1, [20 1]);
[~, startms] = max(h);  if startms < 2, startms = startms + 20; end
nav1 = startms - 1;
nSym = floor((numel(Pi) - startms + 1)/20);
symSoft = sum(reshape(Pi(startms:startms+20*nSym-1), 20, nSym), 1);
sym = symSoft > 0;                               % 1 -> logic 1

sync = [1 1 1 0 1 0 1 1 1 0 0 1 0 0 0 0];
sb = 2*sym - 1;  sp = 2*sync - 1;
corr = conv(sb, fliplr(sp), 'valid');            % corr(k) = sum sb(k:k+15).*sp
cand = find(abs(corr) == 16);
trellis = poly2trellis(7, [171 133]);
crcDet  = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);
gpsPi = 3.1415926535898;  gammaS = (2492.028e6/1176.45e6)^2;
got1 = false; got2 = false; gotT = false; first = NaN;
for k = cand
    if k + 599 > nSym, break; end
    if ~(k + 600 + 15 <= nSym && abs(corr(k+600)) == 16) && ~(k - 600 >= 1 && abs(corr(k-600)) == 16)
        continue                                 % not confirmed by a neighbouring sync word
    end
    blk = sym(k:k+599);
    if corr(k) < 0, blk = ~blk; end              % inverted (Costas 180 deg ambiguity)
    dein = matdeintrlv(double(blk(17:600)), 73, 8);
    dec  = vitdec(dein, trellis, 35, 'term', 'hard');
    [~, err] = crcDet(dec(1:286).');  reset(crcDet);
    if err, continue; end
    b = dec(1:262);
    towc = u_(b(9:25));  sfid = u_(b(28:29)) + 1;
    tStart = 12*towc - 12;
    if isnan(first)
        first = k;  e.TOW1 = tStart;  sfb1 = k - 1;
    end
    switch sfid
        case 1
            e.weeknum = u_(b(31:40));            % IRNSS week (mod 1024)
            e.af0 = s_(b(41:62))*2^-31;  e.af1 = s_(b(63:78))*2^-43;  e.af2 = s_(b(79:86))*2^-55;
            e.N   = u_(b(87:90));        e.toc = u_(b(91:106))*16;
            e.TGD = gammaS*s_(b(107:114))*2^-31;
            e.deltan = s_(b(115:136))*2^-41*gpsPi;
            e.IODC = u_(b(137:144));
            e.health = u_(b(155));                % L5 flag
            e.Cuc = s_(b(157:171))*2^-28;  e.Cus = s_(b(172:186))*2^-28;
            e.Cic = s_(b(187:201))*2^-28;  e.Cis = s_(b(202:216))*2^-28;
            e.Crc = s_(b(217:231))*2^-4;   e.Crs = s_(b(232:246))*2^-4;
            e.idot = s_(b(247:260))*2^-43*gpsPi;
            got1 = true;
        case 2
            e.M0 = s_(b(31:62))*2^-31*gpsPi;   e.toe = u_(b(63:78))*16;
            e.ecc = u_(b(79:110))*2^-33;       e.sqrta = u_(b(111:142))*2^-19;
            e.omegae = s_(b(143:174))*2^-31*gpsPi;  e.w = s_(b(175:206))*2^-31*gpsPi;
            e.omegadot = s_(b(207:228))*2^-41*gpsPi;  e.i0 = s_(b(229:260))*2^-31*gpsPi;
            got2 = true;
        otherwise
            if u_(b(31:36)) == 9 && u_(b(191:193)) == 0     % MT 9, GNSS ID 0 = GPS
                e.igtA0 = s_(b(129:144))*2^-35;  e.igtA1 = s_(b(145:157))*2^-51;
                e.igtA2 = s_(b(158:164))*2^-68;  e.igtTot = u_(b(165:180))*16;
                e.igtWNot = u_(b(181:190));
                gotT = true;
            elseif u_(b(31:36)) == 11 && ~isfield(e,'navIonoAlpha')   % MT 11: ionospheric coefficients
                q = @(a) s_(b(a:a+7));
                e.navIonoAlpha = [q(175)*2^-30, q(183)*2^-27/gpsPi, q(191)*2^-24/gpsPi^2, q(199)*2^-24/gpsPi^3];
                e.navIonoBeta  = [q(207)*2^11,  q(215)*2^14/gpsPi,  q(223)*2^16/gpsPi^2,  q(231)*2^16/gpsPi^3];
            end
    end
    if got1 && got2 && gotT && isfield(e,'navIonoAlpha'), break; end
end
if got1 && got2
    e.IODE2 = e.IODC;  e.IODE3 = e.IODC;
    e.updateflag = double(e.health == 0);
    e.igtValid = gotT;
    if ~gotT     % no MT 9 from this satellite: SDR_main takes it from another NavIC channel
        e.igtA0 = 0; e.igtA1 = 0; e.igtA2 = 0; e.igtTot = e.toe; e.igtWNot = e.weeknum;
    end
    fprintf('NavIC PRN %2d: SF1/SF2 decoded (IODEC %d, toe %d), MT9 %d: IRNWT-GPST %.2f ns, TOW %d at symbol %d\n', ...
        prn, e.IODC, e.toe, gotT, e.igtA0*1e9, e.TOW1, sfb1);
else
    fprintf('NavIC PRN %2d: ephemeris incomplete (SF1 %d SF2 %d)\n', prn, got1, got2);
end
end

function v = u_(bits)
v = sum(double(bits(:)') .* 2.^(numel(bits)-1:-1:0));
end

function v = s_(bits)
n = numel(bits);  v = u_(bits);
if bits(1), v = v - 2^n; end
end
