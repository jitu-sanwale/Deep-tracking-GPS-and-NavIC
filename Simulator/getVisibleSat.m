function [satList,elevation,azimuth] = getVisibleSat(eph,startTime,RxPos,settings)
% -------------------------------------------------------------------------
%                  SoftSim: GPS IF signal simulator 
% Author: 
%        Yafeng Li 
%        @ Beijing Information Science and Technology University(BISTU)
% 2021. 02. 18
% -------------------------------------------------------------------------
%
nEph = numel(eph);                % 32 GPS PRNs or 14 NavIC PRNs
elevation = nan(1,nEph);
azimuth = nan(1,nEph);

for svIndex = 1:nEph
    if isempty(eph(svIndex).sqrta) || eph(svIndex).sqrta == 0
        continue                       % PRN without a broadcast ephemeris
    end
    satPos = getSatPos(startTime,eph(svIndex));
    %--- Calculate distance -----------------------------------------
    rho = (satPos(1) - RxPos(1))^2 + (satPos(2) - RxPos(2))^2 + ...
        (satPos(3) - RxPos(3))^2;
    traveltime = sqrt(rho) / settings.c ;
    
    %--- Correct satellite position (do to earth rotation) --------
    % Convert SV position at signal transmitting time to position
    % at signal receiving time. ECEF always changes with time as
    % earth rotates.
    Rot_X = e_r_corr(traveltime, satPos);
    
    %--- Find the elevation angel of the satellite ----------------
    [azimuth(svIndex), elevation(svIndex), ~] = topocent(RxPos, Rot_X - RxPos);
end

satList = find(elevation>settings.elevationMask);
azimuth = azimuth(satList);
elevation = elevation(satList);


