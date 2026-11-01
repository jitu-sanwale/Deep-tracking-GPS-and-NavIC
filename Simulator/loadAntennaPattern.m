function pattern = loadAntennaPattern(csvFile)
%LOADANTENNAPATTERN Load a GNSS antenna gain pattern from CSV.
%   Expects columns [Azimuth(deg), Elevation(deg), Gain_dBi] on a regular
%   grid, sorted with azimuth varying slowest (as exported by antenna
%   pattern measurement/simulation tools).
%
%   pattern = loadAntennaPattern(csvFile)
%
%   Output fields:
%       azGrid, elGrid  - sorted unique grid axes (deg)
%       gainGrid        - gain (dBi), size [numel(elGrid) x numel(azGrid)]

raw = readmatrix(csvFile);
azCol   = raw(:,1);
elCol   = raw(:,2);
gainCol = raw(:,3);

pattern.azGrid   = unique(azCol);
pattern.elGrid   = unique(elCol);
pattern.gainGrid = reshape(gainCol, numel(pattern.elGrid), numel(pattern.azGrid));
