function [symbols, refTime] = navic_generateNavMsg(payloadTemplate, second)
%NAVIC_GENERATENAVMSG Turn a static 4x292 NavIC subframe payload template
%(from navic_eph2sbf.m) into one 48 s master frame of transmitted symbols.
%
%   [symbols, refTime] = navic_generateNavMsg(payloadTemplate, second)
%
%   Inputs:
%       payloadTemplate - 4x292 matrix from navic_eph2sbf.m (bits 9:25 and
%                         263:292 are placeholders, overwritten here)
%       second          - current simulation time [[s], seconds of week]
%   Outputs:
%       symbols  - 4x600 matrix of bipolar (+-1) channel symbols, one row
%                  per subframe (SF1..SF4), ready to be repeated 20x each
%                  (20 ms/symbol @ 50 sps) onto the carrier.
%       refTime  - start time [s] of the 48 s master frame that was
%                  generated (multiple of 48); TOWC in subframe k encodes
%                  the start of the *next* subframe, refTime + 12*k.
%
% Per IRNSS SPS ICD v1.1: CRC-24Q over bits 1:262 -> bits 263:286, 6 zero
% tail bits (287:292), rate-1/2 K=7 (171/133 octal) convolutional code,
% 73x8 block interleaver, 16-bit EB90 sync word prepended (uncoded). This
% mirrors NavIC_Software_Receiver/include/ephemeris.m and NAVdecoding.m
% exactly (same trellis, same CRC polynomial, same interleaver dimensions).

persistent trellis crcGen syncBits
if isempty(trellis)
    trellis = poly2trellis(7, [171 133]);
    crcGen = comm.CRCGenerator([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);
    syncBits = [1 1 1 0 1 0 1 1 1 0 0 1 0 0 0 0]; % 0xEB90, MSB first
end

refTime = floor(second / 48) * 48;

symbols = zeros(4, 600);
for k = 1:4
    bits292 = payloadTemplate(k, :);

    % TOWC (bits 9:25): count of 12s intervals at which the NEXT subframe
    % starts, i.e. (refTime + 12*k)/12.
    towc = refTime / 12 + k;
    bits292(9:25) = ubits17_(towc);

    % CRC-24Q over the 262-bit information field -> bits 1:286.
    codeword = crcGen(bits292(1:262).');
    bits292(1:286) = codeword.';
    bits292(287:292) = 0; % tail bits

    encoded = convenc(bits292, trellis);          % 292 -> 584
    interleaved = matintrlv(encoded, 73, 8);        % 73x8 block interleave

    symbols(k, :) = [syncBits, interleaved] * 2 - 1; % 0/1 -> bipolar
end

reset(crcGen);

end

function bits = ubits17_(u)
u = mod(round(u), 2^17);
bits = double(bitget(uint32(u), 17:-1:1));
end
