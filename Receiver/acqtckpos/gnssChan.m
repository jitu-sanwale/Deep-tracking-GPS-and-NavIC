function ch = gnssChan(idx, signal, file)
%GNSSCHAN Per-channel signal description for the GPS + NavIC receiver.
%   ch = gnssChan(idx, signal, file)
%   Channels share one index space: idx = 1..32 is GPS L1 C/A PRN idx,
%   idx = 101..114 is NavIC L5 SPS PRN idx-100 (cmn/navic offset 100), so
%   the per-PRN arrays of the receiver (TckResult(idx), eph(idx), sbf.*(idx))
%   hold both systems. Fields:
%     isNavIC  - true for a NavIC channel
%     prn      - PRN within its own system
%     code     - 1x1023 ranging code (+/-1)
%     fc       - carrier frequency [Hz]; lam - carrier wavelength [m]
%     fid      - file id of the IF file of the channel's band
%     ionoScale- factor applied to the L1 Klobuchar delay ((fL1/fc)^2)
%     sys      - 1 GPS, 2 NavIC (column of the inter-system bias)
c = 299792458;
ch.isNavIC = idx > 100;
if ch.isNavIC
    ch.prn  = idx - 100;
    ch.code = generateCAcode_L5_band(ch.prn);
    ch.fc   = 1176.45e6;
    if nargin >= 3, ch.fid = file.fidN; end
    ch.sys  = 2;
else
    ch.prn  = idx;
    ch.code = generateCAcode(idx);
    ch.fc   = signal.Fc;
    if nargin >= 3, ch.fid = file.fid; end
    ch.sys  = 1;
end
ch.lam       = c/ch.fc;
ch.ionoScale = (1575.42e6/ch.fc)^2;
end
