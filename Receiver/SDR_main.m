%Purpose:
%   Main function of the GPS software-defined receiver (SDR) based on 
%   vector tracking
%
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU
% Modified for 8-bit complex IQ, zero-IF, 5 MHz data and continuous
% processing of the whole data file (CT and VT), with result plots.


%%
clc;
clear;
format long g;
codeDir = fileparts(mfilename('fullpath'));
if isempty(codeDir), codeDir = pwd; end
addpath(fullfile(codeDir,'geo'))        % Geo-related functions, e.g. ionospheric correction function
addpath(fullfile(codeDir,'acqtckpos'))  % Acquisition, tracking, and postiong calculation functions
addpath(fullfile(codeDir,'common'))     % trajectory spline, antenna pattern, profile list (shared with the simulator)
addpath(codeDir)
% caches and results go to Receiver\output (or the folder in VT_OUTDIR)
runDir = getenv('VT_OUTDIR'); if isempty(runDir), runDir = fullfile(codeDir, 'output'); end
if ~exist(runDir, 'dir'), mkdir(runDir); end
cd(runDir);
tStart = tic;


%% Parameter initialization 
[file, signal, acq, track, solu, cmn] = initParameters();
rxModes = rxModeSetting(cmn.rxModeDefault);   % [classical vector ultratight]
if file.dataType == 2, dataKind = 'complex I/Q'; else, dataKind = 'real'; end
fprintf('Data file: %s\n  %d ms of data, Fs = %.3f MHz, IF = %.3f MHz, %s\n', file.fileRoute, ...
        file.totalMs, signal.Fs/1e6, signal.IF/1e6, dataKind);


%% Acquisition 
if ~exist(['Acquired_',file.fileName,'.mat'],'file')
    Acquired = acquisition(file,signal,acq);
    save(['Acquired_',file.fileName,'.mat'],'Acquired');    
else
    load(['Acquired_',file.fileName,'.mat']);
end
fprintf('Acquisition Completed. \n');

 
%% Do conventional signal tracking and obtain satellites ephemeris
  if  ~exist(['eph_',file.fileName,'.mat'],'file') || ...
          ~exist(['TckResultCT_forEph_',file.fileName,'.mat'],'file') || ...
            ~exist(['sbf_',file.fileName,'.mat'],'file')
    
    % tracking using conventional DLL and PLL
    TckResultCT_forEph = trackingCT(file,signal,track,Acquired);
    
    % navigaion data decode
    [eph, TckResultCT_forEph, sbf] = naviDecode(Acquired, TckResultCT_forEph);
    
    save(['eph_',file.fileName,'.mat'],'eph'); % ephemeris
    save(['sbf_',file.fileName,'.mat'],'sbf'); % first navi bit point and the beginning of subframe 1
    save(['TckResultCT_forEph_',file.fileName,'.mat'], 'TckResultCT_forEph');
    
else
    load(['TckResultCT_forEph_',file.fileName,'.mat']);
    load(['eph_',file.fileName,'.mat']);
    load(['sbf_',file.fileName,'.mat']);
end
% sbf files decoded before the naviDecode sfb1 off-by-one fix: correct once
if ~isfield(sbf,'sfb1ZeroBased')
    sbf.sfb1 = sbf.sfb1 - 1;  sbf.sfb1ZeroBased = true;
    save(['sbf_',file.fileName,'.mat'],'sbf');
    fprintf('Corrected cached sbf_%s.mat to the 0-based subframe-1 bit index.\n', file.fileName);
end


%% Find satellites that can be used to calculate user position
 posSV  = findPosSV(file,Acquired,eph);
 load(['nAcquired_',file.fileName,'.mat']);
 Acquired = nAcquired;
 % ephemeris sanity check: a decoded set that passed parity can still be
 % inconsistent. IS-GPS-200: the 8 LSBs of IODC (subframe 1) must equal IODE in
 % subframes 2 and 3 for the three to form one set. catapult_trap: PRN 30 paired
 % a corrupt subframe 1 (IODC 71, af0 off by 220 us = 66 km) with subframes 2/3
 % of IODE 86 - the start fix landed 25 km off and the classical EKF went NaN.
 % Also checks toe/toc within the week and a plausible orbit radius.
 okEph = true(1, numel(Acquired.sv));
 for k = 1:numel(Acquired.sv)
     prn = Acquired.sv(k);
     e = eph(prn);
     okEph(k) = mod(e.IODC(1),256) == e.IODE2(1) && e.IODE2(1) == e.IODE3(1) && ...
                e.toe(1) >= 0 && e.toe(1) < 604800 && e.toc(1) >= 0 && e.toc(1) < 604800;
     if okEph(k)
         rs = svPosVel(prn, eph, e.TOW1(1), 1);
         if prn > 100   % NavIC GSO/IGSO (a ~ 42164 km)
             okEph(k) = all(isfinite(rs)) && norm(rs) > 41.0e6 && norm(rs) < 43.5e6;
         else           % GPS MEO
             okEph(k) = all(isfinite(rs)) && norm(rs) > 25.5e6 && norm(rs) < 27.5e6;
         end
     end
     if ~okEph(k), fprintf('PRN %d: decoded ephemeris fails the sanity check - not used.\n', prn); end
 end
 for f = fieldnames(Acquired)'
     Acquired.(f{1}) = Acquired.(f{1})(okEph);
 end
 if length(Acquired.sv) < 4
     error('Only %d satellites with complete ephemeris - cannot compute a position.', length(Acquired.sv));
 end
 fprintf('Satellites used for positioning: %s\n', mat2str(Acquired.sv));
 nGPS = sum(Acquired.sv < 100);
 if cmn.useGPS && nGPS < 4
     error('Only %d GPS satellites - the receiver clock is referred to GPS time and needs >= 4.', nGPS);
 end

 % NavIC channels: every pseudorange is formed from TOW1 - (sfb1*20 ms + nav1 +
 % code delay) with one common time origin (pr_est.m). The NavIC subframes
 % start at multiples of 12 s and the GPS subframe 1 at multiples of 30 s,
 % so the NavIC subframe count is re-based to the GPS reference TOW1:
 %   sfb1' = sfb1 + (TOW1_GPS - TOW1_NavIC)/20 ms,  TOW1' = TOW1_GPS,
 % which leaves every transmit time unchanged (both are system times of
 % week; IRNSS time = GPS time + dt_IG of message type 9, applied in svPosVel).
 refPrn = Acquired.sv(find(Acquired.sv < 100, 1));
 if isempty(refPrn), refPrn = Acquired.sv(1); end   % NavIC-only: first NavIC channel
 for prn = Acquired.sv(Acquired.sv > 100)
     dTow = eph(refPrn).TOW1(1) - eph(prn).TOW1(1);
     sbf.sfb1(prn) = sbf.sfb1(prn) + round(dTow/0.02);
     eph(prn).TOW1 = eph(refPrn).TOW1(1);
 end
 % message type 9 (IRNSS-GPS time offset) is a system parameter: a NavIC
 % channel that did not decode it uses the copy of another NavIC channel
 nv = Acquired.sv(Acquired.sv > 100);
 % (ephemeris caches decoded before the igtValid flag: a zero A0 = not decoded)
 if ~isempty(nv) && isfield(eph, 'igtA0')
     if isfield(eph, 'igtValid')
         valid = arrayfun(@(p) ~isempty(eph(p).igtValid) && eph(p).igtValid, nv);
     else
         valid = arrayfun(@(p) ~isempty(eph(p).igtA0) && eph(p).igtA0 ~= 0, nv);
     end
     okT = nv(valid);
     for prn = nv(~valid)
         if ~isempty(okT)
             for f = {'igtA0','igtA1','igtA2','igtTot','igtWNot'}
                 eph(prn).(f{1}) = eph(okT(1)).(f{1});
             end
             fprintf('NavIC PRN %d: message type 9 taken from PRN %d\n', prn-100, okT(1)-100);
         end
     end
 end
 if any(Acquired.sv > 100)
     fprintf('NavIC channels (PRN %s) re-based to the GPS TOW %d s\n', mat2str(Acquired.sv(Acquired.sv>100)-100), eph(refPrn).TOW1(1));
 end


%% Broadcast ionosphere parameters and day of year from the navigation message
global ALPHA BETA
for prn = Acquired.sv(Acquired.sv < 100)      % GPS broadcast Klobuchar (L1-referred)
    if isfield(eph(prn),'ionoAlpha') && ~isempty(eph(prn).ionoAlpha)
        ALPHA = eph(prn).ionoAlpha;  BETA = eph(prn).ionoBeta;
        fprintf('Using broadcast Klobuchar parameters from PRN %d\n', prn);
        break;
    end
end
if all(ALPHA == 0)   % no GPS: Klobuchar coefficients of NavIC message type 11 (L1-referred here)
    for prn = Acquired.sv(Acquired.sv > 100)
        if isfield(eph(prn),'navIonoAlpha') && ~isempty(eph(prn).navIonoAlpha)
            ALPHA = eph(prn).navIonoAlpha;  BETA = eph(prn).navIonoBeta;
            fprintf('Using NavIC message type 11 ionospheric coefficients from PRN %d\n', prn-100);
            break;
        end
    end
end
wn10     = mod(eph(Acquired.sv(1)).weeknum(1), 1024);
gpsWeek  = wn10 + 2048;                    % assumes data after the April 2019 week rollover
gpsDays  = gpsWeek*7 + eph(Acquired.sv(1)).TOW1(1)/86400;
dn       = datenum(1980,1,6) + gpsDays;
dv       = datevec(dn);
cmn.doy  = floor(dn - datenum(dv(1),1,1)) + 1;
fprintf('GPS week %d, TOW %.0f s -> %s (day of year %d)\n', gpsWeek, eph(Acquired.sv(1)).TOW1(1), datestr(dn), cmn.doy);


%% Positioning over the whole file: classical (always - VT and UT are handed
%  over from it), then vector and/or ultratight tracking
[TckResultCT, navSolutionsCT] = trackingCT_POS(file,signal,track,cmn,Acquired, solu.cnslxyz, eph, sbf);
load (['tckRstCT_',file.fileName,'.mat']);   % CN0_CT

if rxModes(2)
    try
        % VT = the ultratight VDLL + CaP framework without INS aiding
        [TckResultVT, navSolutionsVT] = trackingUT_POS(file,signal,track,cmn,Acquired,navSolutionsCT.refPos,eph,sbf,TckResultCT,navSolutionsCT,[],'none');
        load (['tckRstVT_',file.fileName,'.mat']);   % CN0_VT
    catch errVT
        % a diverged VT must not stop the ultratight run
        fprintf('Vector tracking aborted: %s\n', errVT.message);
        navSolutionsVT = []; CN0_VT = []; TckResultVT = [];
    end
else
    navSolutionsVT = []; CN0_VT = []; TckResultVT = [];
end

imu = [];
if rxModes(3)
    if ~exist(cmn.imuFile,'file')
        error('IMU file not found: %s (generate the profile with the simulator)', cmn.imuFile);
    end
    load(cmn.imuFile, 'imu');
    fprintf('IMU: %s, %d samples @ %g Hz, accel range +/-%g g\n', cmn.imuFile, numel(imu.t), imu.fs, imu.accelRangeG);
    [TckResultUT, navSolutionsUT] = trackingUT_POS(file,signal,track,cmn,Acquired,navSolutionsCT.refPos,eph,sbf,TckResultCT,navSolutionsCT,imu);
else
    navSolutionsUT = []; TckResultUT = [];
end

fprintf('Tracking and Positioing Completed in %.1f min.\n\n', toc(tStart)/60);


%% Plot and save results
outDir = fullfile(runDir, ['results_', file.fileName]);
plotNavResults(navSolutionsCT, navSolutionsVT, CN0_CT, CN0_VT, Acquired.sv, file, outDir);
if exist(cmn.truthFile,'file')
    navCTcmp = navSolutionsCT; if ~rxModes(1), navCTcmp = []; end
    truthSummary = compareTruth(navCTcmp, navSolutionsVT, navSolutionsUT, file, cmn, outDir, imu);
end
fclose(file.fid);
if isfield(file,'fidN') && file.fidN > 0, fclose(file.fidN); end
