function gainDb = lookupAntennaGain(pattern, azDeg, elDeg, belowSlopeDbPerDeg, belowFloorDb)
%LOOKUPANTENNAGAIN Interpolate antenna gain (dBi) from a pattern loaded
%   by loadAntennaPattern.m.
%
%   gainDb = lookupAntennaGain(pattern, azDeg, elDeg)
%   gainDb = lookupAntennaGain(pattern, azDeg, elDeg, belowSlopeDbPerDeg, belowFloorDb)
%
%   azDeg is wrapped to the pattern's azimuth range and elDeg is clamped
%   to the pattern's elevation/off-boresight range before lookup, so
%   out-of-table queries degrade to the nearest edge instead of NaN.
%
%   With the optional roll-off arguments, directions BELOW the table's
%   lowest elevation (for the receive antenna: below the antenna's own
%   horizon, e.g. satellites "above" an inverted aircraft) are attenuated
%   by belowSlopeDbPerDeg per degree below that edge, down to belowFloorDb
%   relative to the edge gain - a generic patch back-lobe model. Without
%   them the edge value is held, which made an inverted antenna lose no
%   signal at all.

az = mod(azDeg, 360);
el = min(max(elDeg, pattern.elGrid(1)), pattern.elGrid(end));
gainDb = interp2(pattern.azGrid, pattern.elGrid, pattern.gainGrid, az, el, 'linear');
if nargin >= 5 && elDeg < pattern.elGrid(1)
    gainDb = gainDb + max(-belowSlopeDbPerDeg*(pattern.elGrid(1) - elDeg), belowFloorDb);
end
