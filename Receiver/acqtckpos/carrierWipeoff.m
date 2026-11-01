function [I, Q] = carrierWipeoff(raw, Wave, signal)
%Purpose:
%   Carrier wipe-off. Wave is the replica phase (rad).
%   mixSign = -1 : multiply by exp(-j*Wave) -> correct for complex baseband/IF
%                  data with a non-inverted spectrum (I = real, Q = imag)
%   otherwise    : original GPSSDR_vt convention (multiply by exp(+j*Wave),
%                  I = imag, Q = real), used for real IF data
if isfield(signal,'mixSign') && signal.mixSign < 0
    bb = raw .* exp(-1i.*Wave);
    I  = real(bb);
    Q  = imag(bb);
else
    bb = raw .* exp(1i.*Wave);
    I  = imag(bb);
    Q  = real(bb);
end
end
