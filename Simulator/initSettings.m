function settings = initSettings()
% -------------------------------------------------------------------------
%                   SoftSim: GPS IF signal simulator 
% Author: 
%        Yafeng Li 
%    @ Beijing Information Science and Technology University(BISTU)
%    2022. 08. 18
% -------------------------------------------------------------------------
%
%% Simulator settings ============================================
% Ground-truth receiver trajectory: [time(s), lat(deg), lon(deg), alt(m),
% roll(deg), pitch(deg), heading(deg)]. Prompted at run time by
% signalGenerating.m (via TrajectorySelectionSetting.m), same pattern as
% multipath/filter/chirp/iono-tropo below - the menu lists every *.csv in
% this folder plus an option to (re)generate the high-dynamics
% jerk-doublet stress-test trajectory (HighDynamicsTrajectory.csv, jerk
% values [0 20 40 80 100 120] m/s^3 as successive doublet maneuvers along
% one continuous trajectory - see generateHighDynamicsTrajectory.m). Set
% defaultTrajFile to a filename (e.g. 'static_Trajectory.csv') to skip
% the prompt for scripted runs; leave empty ([]) to always be prompted.
settings.defaultTrajFile      = [];
% Restrict signal generation to a single (highest-elevation) visible
% satellite - intended for the high-dynamics trajectory, to isolate one
% tracking channel's response from other satellites' independently
% varying geometry/dynamics. Prompted alongside trajFile above;
% TrajectorySelectionSetting.m suggests 1 by default when
% HighDynamicsTrajectory.csv is chosen and 0 otherwise. Set
% defaultSingleSatEnable to 1/0 to skip that prompt.
settings.defaultSingleSatEnable = [];
% Signal length(milliseconds) to be simulated. projectile_Vivek.csv spans
% 0-54.3s; kept just inside that range so no extrapolation is needed.
settings.msToProcess        = 60000;%200000;%91000%;54000;%;200000;%54000;        %[ms]
% Ephemeris (RINEX navigation) file selection. Choose 'brdc' for the
% RINEX 2 GPS-only broadcast file, or 'IITK' for the RINEX 3
% mixed-constellation file (rinexeV3.m parses GPS records only; other
% systems present in the file - GLONASS/Galileo/BeiDou/QZSS/IRNSS - are
% skipped). Both are parsed into the same eph/ionoutc format, so nothing
% else needs to change when switching between them.
settings.ephemerisSource = 'IITK';   % 'brdc' or 'IITK'
switch settings.ephemerisSource
    case 'brdc'
        settings.rinexVersion = 2;
        settings.rinexfile    = 'brdc0500.26n';%'brdc0530.26n';%'brdc3540.14n';%
    case 'IITK'
        settings.rinexVersion = 3;
        settings.rinexfile    = 'IITK00IND_R_20250011000_01H_MN.rnx';
    otherwise
        error(['initSettings: unknown settings.ephemerisSource ''%s'' ' ...
            '(expected ''brdc'' or ''IITK'').'], settings.ephemerisSource);
end
% Name of the IF file to store the generated signals
settings.IfFile             = 'gps_0IF_Fs5MHz_8bit_IQ_HighDynamics_JerkDoublet.bin';
% Data type used to store one sample
settings.dataType           = 'int8';
% File Types
%1 - 8 bit real samples S0,S1,S2,...
%2 - 8 bit I/Q samples I0,Q0,I1,Q1,I2,Q2,...
settings.fileType           = 2;
% Intermediate, sampling and code frequencies
settings.IF                 = 0;                    % [Hz]  zero-IF (complex baseband)
settings.samplingFreq       =  5e6;              % [Hz]
settings.codeFreqBasis      = 1.023e6;              % [Hz]
% Define number of chips in a code period
settings.codeLength         = 1023.;
% Nominal carrier frequency
settings.carrFreqBasis      = 1575.42e6;            % [Hz]
% White noise amplitude (Sigma): should be set according to the ADC
% quantization level and max C/No values 
settings.gwnAmp             = 40;                    % 
% Front end filter order
settings.FeOrder            = 20;
% Enable FE filter or not. Prompted at run time by signalGenerating.m
% (via FilterEnableSetting.m), same as multipath/chirp/iono-tropo below.
% Set defaultFilterEnable to 1/0 to skip the prompt for scripted runs;
% leave empty ([]) to always be prompted.
settings.defaultFilterEnable = 0;
% Front end bandwidth [Hz]. Prompted at run time by signalGenerating.m
% (via FrontEndBandwidthSetting.m, only when the filter is enabled),
% which enforces it stays within the range the current
% samplingFreq/IF/codeFreqBasis can support and suggests a default. Set
% defaultFeBandwidth below to skip the prompt.
settings.defaultFeBandwidth = [];
% Elevation mask to exclude signals from satellites at low elevation
settings.elevationMask      = 5;%10;           %[degrees 0 - 90]
% Uniform C/No (dB-Hz) applied to all visible satellites, skipping the
% interactive CNoSetting() prompt so runs can be scripted. Set to [] to
% get the original interactive prompt instead.
settings.defaultCNoDB       = 45;

%% Multipath settings =============================================
% Statistical multipath model from:
%   Wang, Y.; Chen, X.; Liu, P. "Statistical Multipath Model Based on
%   Experimental GNSS Data in Static Urban Canyon Environment."
%   Sensors 2018, 18, 1149. https://doi.org/10.3390/s18041149
% Leave both empty ([]) so signalGenerating.m always prompts at run time
% for whether to enable multipath and how many rays per satellite to
% simulate. Set explicit values here to skip the prompts for scripted
% runs (e.g. defaultMultipathEnable=1, defaultNumMultipathRays=5).
settings.defaultMultipathEnable    = [];
settings.defaultNumMultipathRays   = [];

%% Antenna gain pattern settings ===================================
% Rx: receiver antenna gain vs. satellite elevation/azimuth as seen from
% the receiver's local horizon. Tx: satellite transmit antenna gain vs.
% the receiver's off-boresight (nadir) angle as seen from the satellite
% (boresight assumed Earth-center-pointing). Both files are
% [Azimuth(deg), Elevation(deg), Gain_dBi] grids.
settings.rxAntennaEnable    = 1;            % 1-on, 0-off
settings.rxAntennaFile      = 'Rx_GNSS_Antenna_Pattern.csv';
% Rx gain for satellites BELOW the antenna's own horizon (the pattern file
% stops at 0 deg): roll-off per degree below it, down to a floor relative
% to the 0-deg gain (generic patch back lobe; matters once the platform
% banks steeply or flies inverted)
settings.rxBelowHorizonSlopeDb = 0.5;       % dB/deg
settings.rxBelowHorizonFloorDb = -20;       % dB
settings.txAntennaEnable    = 1;            % 1-on, 0-off
settings.txAntennaFile      = 'Tx_GNSS_Antenna_Pattern.csv';

%% Chirp (sweep) jamming interference ==============================
% Linear-chirp jammer model from:
%   Sun, K.; Elhajj, M.; Ochieng, W.Y. "A GNSS Anti-Interference Method
%   Based on Fractional Fourier Transform." IEEE Trans. Aerosp. Electron.
%   Syst. 2024, 60(5), 5636-5651. doi:10.1109/TAES.2024.3366147
%
% JNR and the sweep frequency deviation are fixed here (not prompted) -
% edit these values directly to change them. Deviation must stay within
% +/- samplingFreq/2 (Nyquist, here +/-2.5 MHz at the default 5 MHz
% sampling rate) or ChirpJammerSetting.m errors at run time.
settings.chirpJnrDB          = 10;      %[dB]  paper tested -4 to 20 dB
settings.chirpDeltaF0        = 2e6;     %[Hz]  left sweep frequency deviation from L1
settings.chirpDeltaF1        = -2e6;    %[Hz]  right sweep frequency deviation from L1
%
% Leave the rest empty ([]) so signalGenerating.m always prompts at run
% time for whether to enable the jammer and its remaining (sweep-period
% and injection-window) parameters. Set explicit values here to skip
% those prompts for scripted runs (e.g. defaultChirpEnable=1,
% defaultChirpSweepPeriodMs=0.5).
% defaultChirpInjStartSec/defaultChirpInjEndSec confine the jammer to an
% elapsed-time window (seconds from the start of the simulation run)
% instead of running for the whole pass; leave empty to be prompted
% (defaults to the full run duration if you just press Enter there).
settings.defaultChirpEnable         = [];
settings.defaultChirpSweepPeriodMs  = [];
settings.defaultChirpInjStartSec    = [];
settings.defaultChirpInjEndSec      = [];

%% Ionospheric and tropospheric delay ===============================
% Ionosphere: Klobuchar model (klobuchar.m) using alpha/beta coefficients
% broadcast in settings.rinexfile's header - the same coefficients also
% get embedded in the generated navigation message (subframe 4 page 18,
% see eph2sbf.m), so a receiver decoding the broadcast message recovers
% the same values used to inject this delay.
% Troposphere: Saastamoinen zenith-delay model (saastamoinenTropo.m)
% with a simple cosecant elevation mapping and standard-atmosphere
% surface weather (deliberately NOT the same model as the receiver's own
% tropospheric correction - Common/tropo.m, Goad & Goodman - so some
% residual error after correction is expected, same as with real
% receivers whose tropo model never exactly matches the true atmosphere).
settings.ionoEnable  = 1;           % 1-on, 0-off
settings.tropoEnable = 1;           % 1-on, 0-off
settings.tropoP0     = 1013.25;     % surface pressure [mbar/hPa]
settings.tropoT0     = 288.15;      % surface temperature [K] (ICAO standard atmosphere, 15 C)
settings.tropoRH     = 50;          % relative humidity [%]

%% Constants ======================================================
settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = 68.802;       %[ms] Initial signal travel time


