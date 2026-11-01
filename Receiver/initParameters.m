function [file, signal, acq, track, solu, cmn] = initParameters()
%Purpose:
%   Parameter initialization
%Inputs: 
%	None
%Outputs:
%	file        - parameters related to the data file to be processed,a structure
%	signal      - parameters related to signals,a structure
%	acq         - parameters related to signal acquisition,a structure
%	track       - parameters related to signal tracking,a structure
%	solu        - parameters related to navigation solution,a structure
%	cmn         - parameters commmonly used,a structure
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU
% Modified for 8-bit complex IQ, zero-IF, 5 MHz data (gps_0IF_Fs5MHz_8bit_IQ_*.bin)


%% File parameters
codeDir             = fileparts(mfilename('fullpath'));
% data folder: holds IFdata\ (IF files) and FighterProfiles\ (truth and IMU),
% the layout written by the simulator. Default: the simulator folder of this
% package; set the environment variable VT_DATADIR to use any other folder.
dataDir             = getenv('VT_DATADIR');
if isempty(dataDir), dataDir = fullfile(codeDir, '..', 'Simulator'); end
file.fileName       = 'gps_0IF_Fs5MHz_8bit_IQ_Fighter_baseline_cruise';
file.fileRoute      = fullfile(dataDir, 'IFdata', [file.fileName,'.bin']);
if ~isempty(getenv('VT_DATAFILE'))          % optional override (used for testing)
    file.fileRoute = getenv('VT_DATAFILE');
    [~, file.fileName] = fileparts(file.fileRoute);
end
file.skip           = 100;  % ms skipped at the start of the file
file.fid            = fopen(file.fileRoute,'r','ieee-le');
if file.fid < 0, error('Cannot open data file: %s', file.fileRoute); end
file.skiptimeVT     = 2000; % when vector tracking begins (ms after start of positioning)
% Constellations: 'GPS' (L1 C/A only), 'NavIC' is never used alone here, or
% 'GPS+NavIC' (default when the NavIC L5 file of the data set exists). The
% NavIC L5 SPS file has the same name with 'gps_' replaced by 'navic_' (or is
% given by VT_DATAFILE_NAVIC); it shares the sample grid of the GPS file.
% Override with the environment variable VT_GNSS.
[fDir, fBase] = fileparts(file.fileRoute);
file.fileRouteN = fullfile(fDir, [regexprep(fBase, '^gps_', 'navic_'), '.bin']);
if ~isempty(getenv('VT_DATAFILE_NAVIC')), file.fileRouteN = getenv('VT_DATAFILE_NAVIC'); end
cmn.gnss = getenv('VT_GNSS');
if isempty(cmn.gnss)
    if exist(file.fileRouteN, 'file'), cmn.gnss = 'GPS+NavIC'; else, cmn.gnss = 'GPS'; end
end
cmn.useNavIC = contains(upper(cmn.gnss), 'NAVIC');
cmn.useGPS   = contains(upper(cmn.gnss), 'GPS');   % 'NavIC' alone: NavIC-only receiver
file.fidN = -1;
if cmn.useNavIC
    file.fidN = fopen(file.fileRouteN, 'r', 'ieee-le');
    if file.fidN < 0, error('Cannot open NavIC data file: %s', file.fileRouteN); end
    if cmn.useGPS
        file.fileName = [file.fileName '_GPS_NavIC'];   % separate caches and results folder
    else
        file.fileName = [file.fileName '_NavIC'];
    end
end
file.dataType       = 2;    % 1:I (real); 2:I/Q (complex, interleaved I,Q)
file.dataPrecision  = 1;    % 1:int8; 2:int16 (bytes per component)


%% Signal parameters
signal.IF               = 0;        % unit: Hz (zero-IF complex baseband)
signal.Fs               = 5e6;      % unit: Hz
signal.Fc               = 1575.42e6; 	
signal.codeFreqBasis	= 1.023e6;  	
signal.ms               = 1e-3; % unit: s
signal.Sample           = ceil(signal.Fs*signal.ms);	
signal.codelength       = signal.codeFreqBasis * signal.ms;
signal.mixSign          = +1;   % spectrum convention of the I/Q data:
                                % -1: I+jQ = s*exp(+j*2*pi*f*t) (standard, e.g. gps-sdr-sim)
                                % +1: inverted spectrum (I-jQ). This file is INVERTED:
                                %     with -1 the measured Dopplers were the exact negatives
                                %     of those predicted from the ephemeris.


%% Data length (derived from file size)
fInfo           = dir(file.fileRoute);
bytesPerMs      = signal.Sample * file.dataType * file.dataPrecision;
file.totalMs    = floor(fInfo.bytes / bytesPerMs);


%% Acquisition parameters
acq.prnList     = 1:32;     % PRN list (GPS L1 C/A)
acq.prnListN    = 1:14;     % NavIC L5 SPS PRNs (channel index 100+PRN, see gnssChan.m)
if ~cmn.useGPS, acq.prnList = []; end   % NavIC-only receiver
acq.freqStep    = 500;      % unit: Hz
acq.freqMin     = -10000;   % Minimum Doppler frequency
acq.freqNum     = 2*abs(acq.freqMin)/acq.freqStep+1;    % number of frequency bins
acq.L           = 10;       % number of ms to perform FFT
acq.threshold   = 20;       % acquisition threshold (dB); 18 gave false detections on 5 MHz IQ test data


%% Tracking parameters
track.CorrelatorSpacing  	= 0.5;  % unit: chip
track.DLLBW               	= 2;	% unit: Hz
track.DLLDamp           	= 0.707; 
track.DLLGain            	= 0.1;	
track.PLLBW              	= 20;  	
track.PLLDamp             	= 0.707;
track.PLLGain              	= 0.25; 	
track.msEph                 = min(50000, file.totalMs - file.skip - 100); % at least ~36 s needed for a full frame
track.msToProcessCT       	= file.totalMs - file.skip - 200;            % whole file (unit: ms)
track.msToProcessVT         = track.msToProcessCT - file.skiptimeVT - 100;
track.pdi                   = 1;
track.progressEvery         = 1000; % print progress every N ms
track.cn0Mask               = 25;   % dB-Hz, channels below are excluded from the EKF update
track.cn0Floor              = 0;    % dB-Hz, C/N0 reported for a lost-signal window (mean(Z)^2 <= var(Z))
track.elMask                = 5;    % deg, elevation mask for the EKF update
track.prGate                = 500;  % m, CT pseudorange innovation gate
track.gateStart             = 2000; % ms, innovation gate active after this epoch

% Ultratight VDCaP (trackingUT_POS.m)
track.utUpdateMs            = 10;   % ms, EKF update interval (prefilter length L)
track.utPLLBW               = 5;    % Hz, CaP loop noise bandwidth (dynamics removed by the INS)
track.utPLIreset            = 0.9;  % CaP loop filter reset when PLI < this after an update
track.utPLIcohMs            = 5;    % ms, coherent sum length for the PLI (divides utUpdateMs)
track.utCodeSteerTau        = 0.2;  % s, VDLL replica steering time constant
track.utGate                = 5;    % normalised innovation gate (sigma)
track.utPrFloor             = 1.0;  % m, pseudorange noise floor
track.utPrrSigma            = 0.05; % m/s, pseudorange-rate noise at 45 dB-Hz
track.utClkPhasePSD         = 1;    % m^2/s, receiver clock phase PSD
track.utClkFreqPSD          = 1;    % m^2/s^3, receiver clock frequency PSD

% Vector tracking VDLL + CaP without INS (trackingUT_POS.m, aiding 'none'):
% the ultratight channels, measurements and filter update, with the INS
% replaced by a Singer position-velocity-acceleration model (Singer, 1970).
% sigma_a is limited by closed-loop stability: the EKF corrects every 10 ms
% from the residual CaP frequency, and on the GPS turn sigma_a = 1 m/s^2 lost
% lock at manoeuvre onset while 10 m/s^2 (and the 48 m/s^2 of Singer's 9 g
% envelope rule) was unstable in level flight; 5 m/s^2 is used for all profiles.
track.vtSingerSigmaA        = 5;    % m/s^2, acceleration standard deviation
track.vtSingerTau           = 1;    % s, manoeuvre (acceleration) correlation time
track.vtP0Acc               = 1;    % m/s^2, initial acceleration std (level flight at handover)
track.vtPLLBW               = 20;   % Hz, VT residual CaP bandwidth (no INS: 5 Hz could not hold lock in level flight)
% After each EKF update, remove the change of the navigation Doppler f_Nav
% from the CaP loop-filter state (carrier NCO continuous across the update)
track.capNavTransfer        = struct('UT', true, 'VT', true);
track.utAttErr0Deg          = [0.01 -0.01 0.02]; % deg, initial roll/pitch/yaw error (datasheet GNSS-aided)
track.utP0AttDeg            = 0.05; % deg, initial attitude sigma
track.utP0Vel               = 1;    % m/s
track.utP0Pos               = 10;   % m
track.utP0Clk               = 10;   % m
track.utP0Drift             = 1;    % m/s
% GPS-NavIC inter-system bias (receiver L5-L1 hardware delay plus the error
% of the broadcast IRNSS-GPS time offset after message type 9 is applied):
% random walk with this PSD, initial sigma isbP0 (all three filters)
track.isbP0                 = 30;   % m
track.isbPSD                = 1e-4; % m^2/s


%% Navigation parameters
% The initial position is now obtained from a least-squares fix at the first
% positioning epoch. Leave iniPos empty to use that fix as the ENU reference,
% or give [lat(rad) lon(rad) h(m)] to use a known reference point.
solu.iniPos  = [];
if isempty(solu.iniPos)
    solu.cnslxyz = [];
else
    solu.cnslxyz = llh2xyz(solu.iniPos);
end
solu.rate  	= 1000; % unit: Hz


%% commonly used parameters
cmn.vtEnable  	= 1;            % 0: disable vector tracking; 1:enable vector tracking (overridden by the mode prompt in SDR_main)
% Receiver mode (prompted in SDR_main via rxModeSetting.m): 1 classical
% (DLL/PLL + EKF), 2 vector (VDLL + EKF), 3 ultratight (VDCaP + 17-state
% EKF + 6-DOF INS), 4 all three (default). Set here or with the environment
% variable VT_RXMODE to skip the prompt; leave empty to be prompted.
cmn.rxModeDefault = [];
% IMU data for the ultratight mode (GPS_Simulator/run_gen_imu_iNAT.m) and the
% truth trajectory for the comparison (spline model of the CSV the IF file was
% generated from; the 50 Hz original gives the identical spline)
cmn.dataDir     = dataDir;
cmn.imuFile     = fullfile(dataDir, 'FighterProfiles', 'Fighter_baseline_cruise_IMU_200Hz.mat');
cmn.truthFile   = fullfile(dataDir, 'FighterProfiles', 'Fighter_baseline_cruise.csv');
% optional overrides (batch runs, e.g. the fighter manoeuvre profiles)
if ~isempty(getenv('VT_IMUFILE')),   cmn.imuFile   = getenv('VT_IMUFILE');   end
if ~isempty(getenv('VT_TRUTHFILE')), cmn.truthFile = getenv('VT_TRUTHFILE'); end
cmn.cSpeed      = 299792458;    % speed of light, [m/s]
cmn.doy         = 1;            % Day of year; overwritten from the decoded GPS week/TOW


%% ionospheric model (overwritten by the broadcast values in subframe 4 page 18 when decoded)
global ALPHA BETA
ALPHA = [0 0 0 0];
BETA  = [0 0 0 0];
