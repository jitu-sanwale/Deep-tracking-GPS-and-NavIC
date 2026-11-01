function signalGenerating(settings)
% Main script generating IF signals
% Input:
%         settings    - IF simulator settings
%      
% -------------------------------------------------------------------------
%                   SoftSim: GPS IF signal simulator 
% Author: 
%        Yafeng Li 
%    @ Beijing Information Science and Technology University(BISTU)
%    2022. 08. 18
% -------------------------------------------------------------------------
%
%
%% Load trajectory ================================================
% Trajectory CSV and single-satellite restriction (for tracking-loop
% dynamic-stress testing) are selected here, same prompt-based pattern as
% multipath/filter/chirp/iono-tropo below, so the choice takes effect for
% the whole run (settings.trajFile is read further down; singleSatEnable
% is applied once the visible-satellite list is known).
[settings.trajFile, settings.singleSatEnable] = TrajectorySelectionSetting( ...
    settings.defaultTrajFile, settings.defaultSingleSatEnable);

% Sample counts within each transmitting time calculation step (1 ms)
blockSize = round(settings.samplingFreq * 0.001);
% Time interval for each transmitting time calculation step
blockTime = blockSize/settings.samplingFreq;

% Ground-truth trajectory CSV, header row expected:
%   Time [sec], Latitude [deg], Longitude [deg], Height [meter],
%   Roll Angle [deg], Pitch Angle [deg], Heading [deg]
% The three attitude columns are optional; if absent, attitude defaults
% to zero (antenna boresight pointing at local zenith - the simulator's
% previous no-attitude behavior).
trajRaw  = readmatrix(settings.trajFile, 'NumHeaderLines', 1);
trajTime = trajRaw(:,1);
trajLat  = trajRaw(:,2);
trajLon  = trajRaw(:,3);
trajAlt  = trajRaw(:,4);
if size(trajRaw,2) >= 7
    trajRoll    = trajRaw(:,5);
    trajPitch   = trajRaw(:,6);
    trajHeading = trajRaw(:,7);
else
    trajRoll    = zeros(size(trajTime));
    trajPitch   = zeros(size(trajTime));
    trajHeading = zeros(size(trajTime));
end

% If the requested run length exceeds the logged trajectory, the receiver
% would otherwise sit at the trajectory's FINAL point for the remainder
% of the run (qTime is clamped to trajTime(end) just below) - i.e. mostly
% static output with no warning. That silent freeze is easy to trigger
% by mismatch (e.g. leaving settings.msToProcess at a previous
% trajectory's length) and easy to misread as "the trajectory itself is
% static" when it is really just msToProcess outliving it. Clamp instead
% of extrapolating, and say so loudly.
if settings.msToProcess/1000 > trajTime(end)
    warning(['signalGenerating: settings.msToProcess (%.1f s) exceeds the ' ...
        'duration of trajectory file %s (%.1f s) - the receiver would sit ' ...
        'static at the trajectory''s final point for the last %.1f s of the ' ...
        'run. Clamping settings.msToProcess to %.1f s to match the trajectory.'], ...
        settings.msToProcess/1000, settings.trajFile, trajTime(end), ...
        settings.msToProcess/1000 - trajTime(end), trajTime(end));
    settings.msToProcess = floor(trajTime(end)*1000);
end
% Iteration count (one loop per 1 ms of simulated signal)
iterCnt = round(settings.msToProcess/1000/blockTime);

% Resample the trajectory onto the simulator's per-block time grid.
% Clamp beyond the last logged sample so a longer msToProcess just holds
% the final trajectory point instead of extrapolating garbage. Heading is
% interpolated linearly rather than circularly, so this can misbehave
% right at a 360/0 deg wraparound between two logged samples.
qTime    = min((0:iterCnt) * blockTime, trajTime(end));
qLat     = interp1(trajTime, trajLat, qTime, 'linear');
qLon     = interp1(trajTime, trajLon, qTime, 'linear');
qAlt     = interp1(trajTime, trajAlt, qTime, 'linear');
qRoll    = interp1(trajTime, trajRoll, qTime, 'linear');
qPitch   = interp1(trajTime, trajPitch, qTime, 'linear');
qHeading = interp1(trajTime, trajHeading, qTime, 'linear');

% Columns: [lat, lon, alt, ECEFx, ECEFy, ECEFz, roll, pitch, heading]
trajectory = zeros(iterCnt + 1, 9);
for index = 1:(iterCnt + 1)
    [ECEFx,ECEFy,ECEFz] = geo2cartd(qLat(index), qLon(index), qAlt(index), 5);
    trajectory(index,1:9) = [qLat(index), qLon(index), qAlt(index), ECEFx, ECEFy, ECEFz, ...
        qRoll(index), qPitch(index), qHeading(index)];
end

%% Read ephemeris file ============================================
if settings.rinexVersion == 2
    [eph,ionoutc] = rinexeV2(settings.rinexfile);
elseif settings.rinexVersion == 3
    [eph,ionoutc] = rinexeV3(settings.rinexfile);
end
% Klobuchar coefficients as broadcast in this RINEX file's header -
% these are also embedded verbatim into subframe 4 page 18 of the
% generated navigation message (eph2sbf.m), so a receiver that decodes
% the broadcast message would recover the same values used here.
ionoAlpha = [ionoutc.alpha0, ionoutc.alpha1, ionoutc.alpha2, ionoutc.alpha3];
ionoBeta  = [ionoutc.beta0,  ionoutc.beta1,  ionoutc.beta2,  ionoutc.beta3];

%% Constellation ===================================================
% settings.gnss = 'GPS' (default, L1 C/A) or 'NavIC' (L5 SPS). A NavIC run
% keeps the GPS time base - the start time is still the GPS reference
% ephemeris' toc - so that a GPS L1 file and a NavIC L5 file generated for
% the same trajectory share one receiver clock and sample grid. NavIC
% system time is aligned with GPS time (IRNSS SPS ICD v1.1, Sec. 5.7).
isNavIC = isfield(settings,'gnss') && strcmpi(settings.gnss,'NavIC');
startTime = eph(1).toc;
gpsWeek   = eph(1).weekNrm;
if isNavIC
    % real IRNSS broadcast ephemeris (I-records of the same RINEX file),
    % per PRN the record with toe closest to the start time
    eph = rinexeIRNSS(settings.rinexfile, startTime);
    % carrier, PRN code and nav message of NavIC L5 SPS; the Klobuchar
    % delay of the GPS coefficients is referred to L5 by (fL1/fL5)^2
    settings.carrFreqBasis = 1176.45e6;
    % NavIC system time (IRNWT) vs GPS time: the signals are stamped in IRNWT,
    %   IRNWT = GPST + dtIG(t),  dtIG = A0 + A1 (t - tot) + A2 (t - tot)^2 + dtErr
    % A0..tot are broadcast in message type 9 (subframe 3); dtErr is the
    % error of that broadcast prediction (not broadcast). rxL5Delay is the
    % receiver's L5-minus-L1 hardware delay, so the receiver sees a GPS/NavIC
    % inter-system bias of c*(rxL5Delay - dtErr) after applying message type 9.
    % The IRNSS week number is the GPS week - 1024 (start epoch 22 Aug 1999).
    igt = struct('A0', 25e-9, 'A1', 2e-14, 'A2', 0, ...
                 'tot', 16*floor((startTime - 3600)/16), 'WNot', gpsWeek - 1024);
    navicDtErr = 3e-9;  navicRxL5Delay = 12e-9;
    if isfield(settings,'navicTime')            % override from the driver
        f = fieldnames(settings.navicTime);
        for k = 1:numel(f)
            switch f{k}
                case 'dtErr',     navicDtErr = settings.navicTime.dtErr;
                case 'rxL5Delay', navicRxL5Delay = settings.navicTime.rxL5Delay;
                otherwise,        igt.(f{k}) = settings.navicTime.(f{k});
            end
        end
    end
    fprintf('NavIC time: IRNWT-GPST = %.1f ns + %.1e s/s (broadcast, MT 9), broadcast error %.1f ns, rx L5-L1 delay %.1f ns\n', ...
        igt.A0*1e9, igt.A1, navicDtErr*1e9, navicRxL5Delay*1e9);
    navicDtIG = @(tg) igt.A0 + igt.A1*check_t(tg - igt.tot) + igt.A2*check_t(tg - igt.tot).^2 + navicDtErr;
    % truth of the time model, for the receiver's truth comparison: the
    % inter-system bias left after message type 9 is c*(rxL5Delay - dtErr)
    navicTime = struct('igt', igt, 'dtErr', navicDtErr, 'rxL5Delay', navicRxL5Delay, ...
        'isbTruth_m', settings.c*(navicRxL5Delay - navicDtErr), 'startTime', startTime);
    [ifDir, ifBase] = fileparts(settings.IfFile);
    save(fullfile(ifDir, [ifBase '_time.mat']), '-struct', 'navicTime');
end
ionoScale = (1575.42e6/settings.carrFreqBasis)^2;   % 1 for GPS L1

%% Determine the visible satellite list ===========================
% Starting time of the simulation, corresponding to the first positions
% in the trajectory (startTime, set above)
% Initail position of the receiver trajectory
RxPosEcef = trajectory(1,4:6);
% Visible satellite list (PRN#)
[satList,elevation,azimuth] = getVisibleSat(eph,startTime,RxPosEcef',settings);

%% Single-satellite restriction ====================================
% For high-dynamics tracking-loop stress testing (see
% TrajectorySelectionSetting.m / generateHighDynamicsTrajectory.m): keep
% only ONE visible satellite, so the generated IF file exercises one
% tracking channel against the injected dynamics without the confound of
% other satellites' independently varying geometry/dynamics.
%
% Default: the single highest-elevation visible satellite (best raw
% CN0/cleanest signal). Set settings.singleSatPRN to a specific PRN
% instead when the test needs a particular line-of-sight geometry - e.g.
% a purely-horizontal (East) trajectory couples almost none of its
% acceleration into a near-zenith satellite's LOS (cos(elevation) -> 0),
% so a low/moderate-elevation satellite whose azimuth is close to the
% direction of travel gives a far larger, more demonstrative LOS-
% projected acceleration. Leave singleSatPRN empty ([]) for the old
% highest-elevation default.
if settings.singleSatEnable
    if isempty(satList)
        error(['signalGenerating: singleSatEnable is set but no satellites ' ...
            'are visible at the trajectory start time.']);
    end
    % Number of satellites to KEEP once the restriction is enabled.
    % Generalizes the original strict single-satellite switch to a chosen
    % COUNT (see SatelliteCountSetting.m): defaults to 1 - the original
    % behavior - whenever settings.numSatellites is unset/empty, so every
    % existing driver script (none of which sets this field) is
    % unaffected. A count >= 2 keeps the N highest-elevation satellites
    % instead of just one, e.g. for a scenario that still wants a real
    % least-squares position fix (needs >= 4) under the same dynamics/C-No
    % stress. It is clamped to however many satellites are actually
    % visible, so a large count (e.g. 99) is a convenient way to mean
    % "don't restrict, use everything visible."
    if isfield(settings, 'numSatellites') && ~isempty(settings.numSatellites)
        nSat = min(settings.numSatellites, numel(satList));
    else
        nSat = 1;
    end

    if nSat == 1 && isfield(settings, 'singleSatPRN') && ~isempty(settings.singleSatPRN)
        bestIdx = find(satList == settings.singleSatPRN, 1);
        if isempty(bestIdx)
            error(['signalGenerating: settings.singleSatPRN = %d is not in the ' ...
                'visible satellite list at the trajectory start time (visible: %s).'], ...
                settings.singleSatPRN, mat2str(satList));
        end
        fprintf('Single-satellite mode: using PRN %d (elevation %.1f deg, azimuth %.1f deg) as requested, out of %d visible.\n', ...
            satList(bestIdx), elevation(bestIdx), azimuth(bestIdx), numel(satList));
    else
        [~, order] = sort(elevation, 'descend');
        bestIdx = order(1:nSat);
        if nSat == 1
            fprintf('Single-satellite mode: restricting to PRN %d (highest elevation, %.1f deg) out of %d visible.\n', ...
                satList(bestIdx), elevation(bestIdx), numel(satList));
        else
            fprintf('Multi-satellite restriction: keeping the %d highest-elevation satellites (PRNs %s) out of %d visible.\n', ...
                nSat, mat2str(satList(bestIdx)'), numel(satList));
        end
    end
    satList   = satList(bestIdx);
    elevation = elevation(bestIdx);
    azimuth   = azimuth(bestIdx);
end

%% Sky plot and C/N0 seting========================================
skyPlot(azimuth',elevation',satList');
% Set C/No: 1.5dB is added for compesating the cross-correlation interference
% between diffierent PRNs
CNoValuesdB = CNoSetting(satList, settings.defaultCNoDB) + 1.5;
% Covert to unite of Hz
CNoValues = 10.^(CNoValuesdB/10);
% Calculate carrier amplitude 
carrAmp = 2 * settings.gwnAmp * sqrt(CNoValues/settings.samplingFreq);
if  settings.fileType == 2
    carrAmp = carrAmp/sqrt(2);
end

%% Multipath settings ==============================================
% Statistical multipath model from Wang, Chen & Liu, "Statistical
% Multipath Model Based on Experimental GNSS Data in Static Urban Canyon
% Environment," Sensors 2018, 18, 1149, doi:10.3390/s18041149.
[settings.multipathEnable, settings.numMultipathRays] = ...
    MultipathSetting(settings.defaultMultipathEnable, settings.defaultNumMultipathRays);

if settings.multipathEnable
    numMP = settings.numMultipathRays;
    % Per-satellite, per-ray parameters, held constant for the whole run
    % (matches the paper's notion of a multipath "event": roughly-fixed
    % delay/attenuation for a section, while the carrier phase difference
    % Phi_k(t) keeps evolving at the drawn fading frequency).
    mpDelaySec = zeros(length(satList), numMP);   % extra propagation delay [s]
    mpAtt      = zeros(length(satList), numMP);   % linear amplitude ratio alpha_k
    mpFadeFreq = zeros(length(satList), numMP);   % fading frequency [Hz]
    mpPhase0   = zeros(length(satList), numMP);   % initial carrier phase [rad]
    for svIndex = 1:length(satList)
        [gammaP, scaleP] = mpDelayParams(elevation(svIndex));  % Table 4, Eq. (7)
        lambdaP = mpFadeLambda(elevation(svIndex));            % Table 8 (MEO), Eq. (11)
        for k = 1:numMP
            delay_m = gamrnd(gammaP, scaleP);
            attDb   = -0.0032*delay_m - 12.3;                  % Eq. (10), Table 7 mean power attenuation
            fadeMag = exprnd(1/lambdaP);
            fadeSign = sign(rand - 0.5);
            if fadeSign == 0
                fadeSign = 1;
            end
            mpDelaySec(svIndex,k) = delay_m / settings.c;
            mpAtt(svIndex,k)      = 10^(attDb/20);
            mpFadeFreq(svIndex,k) = fadeSign * fadeMag;
            mpPhase0(svIndex,k)   = 2*pi*rand;
        end
    end
    fprintf('Multipath enabled: %d ray(s) per satellite (Wang et al., Sensors 2018 statistical model).\n', numMP);
end

%% Front-end filter enable/bandwidth =================================
% Prompted here (rather than only read from the initSettings.m hardcoded
% default) so the choice always takes effect for this run, the same
% pattern every other toggle in this simulator (multipath, chirp jammer,
% iono/tropo) already uses.
settings.filterEn = FilterEnableSetting(settings.defaultFilterEnable);

% Bandwidth is only needed (and only validated) when the filter is
% actually going to be applied - with filterEn=0 the bandwidth is never
% used, so skipping this avoids forcing a value (and its validation
% error()s) that doesn't matter for this run.
if settings.filterEn
    settings.FeBandwidth = FrontEndBandwidthSetting(settings, settings.defaultFeBandwidth);
    fprintf('Front-end filter bandwidth: %.3f MHz.\n', settings.FeBandwidth/1e6);
else
    fprintf('Front-end filter disabled.\n');
end

%% Chirp (sweep) jamming interference ===============================
% Model from Sun, Elhajj & Ochieng, IEEE TAES 2024:
%   j(t) = A_inst(t) * cos(2*pi*integral(finst,0,t) + phi0)
%   finst(t) = f0 + k*t, resetting every sweep period tj
[settings.chirpEnable, settings.chirpJnrDB, settings.chirpSweepPeriodMs, ...
    settings.chirpDeltaF0, settings.chirpDeltaF1, settings.chirpInjStartSec, ...
    settings.chirpInjEndSec] = ChirpJammerSetting(settings, ...
    settings.defaultChirpEnable, settings.defaultChirpSweepPeriodMs, ...
    settings.defaultChirpInjStartSec, settings.defaultChirpInjEndSec);

if settings.chirpEnable
    % Amplitude from JNR. Unlike C/No (a noise-PSD-referenced density,
    % dB-Hz, hence the /samplingFreq in the carrAmp formula above), JNR is
    % a total-power ratio (dB): JNR_lin = P_jammer / P_noise. The WGN
    % noise generated below has power gwnAmp^2 per real/imag component,
    % i.e. total noise power gwnAmp^2 (fileType 1) or 2*gwnAmp^2 (complex,
    % independent I/Q noise). Equating P_jammer to that and solving for
    % the jammer's peak amplitude gives the same closed form for both
    % file types: jamAmp = gwnAmp*sqrt(2*JNR_lin).
    jnrLinear = 10^(settings.chirpJnrDB/10);
    jamAmp = settings.gwnAmp * sqrt(2*jnrLinear);
    % Chirp rate [Hz/s] and a fixed random initial phase for the jammer.
    chirpK = (settings.chirpDeltaF1 - settings.chirpDeltaF0) / (settings.chirpSweepPeriodMs/1000);
    chirpPhi0 = 2*pi*rand - pi;
    fprintf('Chirp jammer enabled: JNR=%.1f dB, sweep period=%.3f ms, deviation=[%.2f, %.2f] MHz, active %.3f-%.3f s.\n', ...
        settings.chirpJnrDB, settings.chirpSweepPeriodMs, settings.chirpDeltaF0/1e6, settings.chirpDeltaF1/1e6, ...
        settings.chirpInjStartSec, settings.chirpInjEndSec);
end

%% Antenna gain patterns ============================================
% Rx pattern is applied vs. the satellite's elevation/azimuth in the
% antenna's own BODY frame - i.e. the local-horizon elevation/azimuth
% (from topocent.m) rotated by the platform's roll/pitch/heading, via
% applyAttitudeToBodyFrame() below, so a tumbling/maneuvering platform's
% antenna gain toward each satellite is modeled correctly instead of
% always assuming a level, zenith-pointing antenna. Tx pattern is applied
% vs. the receiver's off-boresight (nadir) angle as seen from the
% satellite, assuming the satellite antenna boresight points at Earth
% center (standard for GPS) - this is a satellite-attitude effect,
% unrelated to the receiver platform's own attitude, so it is NOT rotated
% by roll/pitch/heading. A proper Tx boresight azimuth would require
% modeling the satellite's own yaw-steering attitude, which is out of
% scope here; the supplied pattern is azimuth-symmetric anyway, so
% azimuth=0 is passed for the Tx lookup.
if settings.rxAntennaEnable
    rxAntPattern = loadAntennaPattern(settings.rxAntennaFile);
end
if settings.txAntennaEnable
    txAntPattern = loadAntennaPattern(settings.txAntennaFile);
end

%% Initialize data bits for each channle ==========================
for svIndex = 1:length(satList)
    PRN = satList(svIndex);
    if isNavIC
        % NavIC L5 SPS message (navic_eph2sbf / navic_generateNavMsg): 48 s
        % master frame of four 600-symbol subframes, CRC-24Q, rate-1/2 K=7
        % FEC, 73x8 interleaver, EB90 sync; regenerated per master frame
        % in the main loop
        navBits(svIndex).template = navic_eph2sbf(eph(PRN), ionoutc, igt);
        navBits(svIndex).refTime  = -Inf;
        navBits(svIndex).dataBit  = [];
        continue
    end
    % Genaerate Nav message according to the ephemeris
    navBits(svIndex).frameMsg = eph2sbf(eph(PRN),ionoutc); %#ok<*SAGROW>
    
    % Time to generate Nav Messages
    navBits(svIndex).ephTime.second = startTime;
    navBits(svIndex).ephTime.weekNrm = eph(1).weekNrm;
    
    % Modulate WN, TOW and CRC into Nav messages
    [navBits(svIndex).dataWord, navBits(svIndex).refTime] = ...
        generateNavMsg(navBits(svIndex).frameMsg,navBits(svIndex).ephTime,1);
    
    % Extract Nav data bits from nav message
    navBits(svIndex).dataBit = zeros(1,30*60);
    for bitindex = 1:60
        navBits(svIndex).dataBit((bitindex-1) * 30 + 1 : bitindex * 30) = ...
            double(bitget(navBits(svIndex).dataWord(bitindex),30:-1:1)) * 2 - 1;
    end
end

%% Generate C/A code table
% Get a vector with the C/A code sampled 1x/chip
for svIndex = 1:length(satList)
    PRN = satList(svIndex);
    if isNavIC
        caCode = generateCAcode_L5_band(PRN);   % ICD Table 7 (+/-1, logic 1 -> +1)
    else
        caCode = generateCAcode(PRN);
    end
    caCodeTable(svIndex,:) = [caCode(end) caCode caCode(1)];
end

%% Allocate variables and spaces ==================================
% blockSize, blockTime and iterCnt were already computed above while
% building the trajectory grid.

% Initialize Rx Time and Rx postions
RxTime(1) = startTime;
RxTime(2) = RxTime(1) + blockTime;
RxPosEcef(1:2,:)  = trajectory(1:2,4:6);
% Local oscilator frequency in rad
localOsFreq = (settings.carrFreqBasis - settings.IF) * 2 * pi;
% RF signal freq in rad
carrFreqRad = settings.carrFreqBasis * 2 * pi;
% Nav bit period
bitPeriod = 0.001 * 20;
% Front end filter coeffficient
% Front end filter coefficient and FIR filter state (carried across ms
% blocks so the streaming filter output is continuous - a fresh filter()
% call per block, with no state passed in, would insert an artificial
% transient at every block boundary). Only needed when filterEn=1 - see
% the front-end filter bandwidth setup above.
if settings.filterEn
    coef = getFliterCoef(settings);
    filterState = zeros(settings.FeOrder,1);
end
% Open the IF file to save senerated samples
[fid, ~] = fopen(settings.IfFile, 'w');

% Start waitbar
hwb = waitbar(0,'IF signal generating ...');
barTimeMs  = round(iterCnt * blockTime * 1000); % [ms]

%% Generate IF siganl =============================================
disp('IF signal generating is undegoing, please wait ...')
for loopCnt =  1:iterCnt
    % wait bar ------------------------------------------------------------
    if (rem(loopCnt, 100) == 0)
        Ln = newline;
        processStatus = ['Generating: ', int2str(loopCnt), ...
            ' ms ', ' of ', int2str(barTimeMs), ' msec'];
        try
            waitbar(loopCnt/barTimeMs,hwb,processStatus);
        catch
            % The progress bar was closed. It is used as a signal
            % to stop, "cancel" processing. Exit.
            disp('Progress bar closed, exiting...');
            return
        end
    end

    % Sum of loacal signals of all visible satellites
    if settings.fileType == 1
        localSigSum = zeros(1,blockSize); % samples of 1ms length
    elseif settings.fileType == 2
        localSigSum = complex(zeros(1,blockSize));
    end

    % Platform attitude for this block (same satellite-independent value
    % used by every channel below).
    rollDeg    = trajectory(loopCnt,7);
    pitchDeg   = trajectory(loopCnt,8);
    headingDeg = trajectory(loopCnt,9);

    % Time-varying C/No profile (optional) ---------------------------------
    % Recomputes carrAmp for THIS block when settings.cnoProfileEnable is
    % set, so C/No can settle/ramp/dip over the run instead of staying
    % fixed at settings.defaultCNoDB for the whole file (see
    % CNoDipSetting.m and the cno31-maneuver / highdyn-cnodip driver
    % scripts). CNoSetting() with a scalar default is a deterministic
    % broadcast (no randomness), so recomputing it every 1ms block is
    % cheap and exactly reproducible.
    %   settings.cnoProfile.shape  - 'linear_hold' (default, used when
    %                                .shape is absent) or 'halfsine_dip'
    %   settings.cnoProfile.tStart - [s] time the departure from
    %                                settings.defaultCNoDB begins
    %   settings.cnoProfile.tEnd   - [s] time the profile reaches its
    %                                final value: for 'linear_hold' this
    %                                is where it reaches lowDB and HOLDS
    %                                there for the rest of the run; for
    %                                'halfsine_dip' this is where it has
    %                                recovered back to defaultCNoDB and
    %                                holds THAT for the rest of the run.
    %   settings.cnoProfile.lowDB  - [dB-Hz] the profile's low point,
    %                                reached at tEnd (linear_hold) or at
    %                                the dip's midpoint (halfsine_dip).
    % Time reference matches the ground-truth trajectory's own clock (t=0
    % at the start of generation), since RxTime(1) at the top of this
    % block equals startTime plus (loopCnt-1)*blockTime.
    if isfield(settings, 'cnoProfileEnable') && settings.cnoProfileEnable
        tNow  = RxTime(1) - startTime;
        prof  = settings.cnoProfile;
        shape = 'linear_hold';
        if isfield(prof, 'shape') && ~isempty(prof.shape)
            shape = prof.shape;
        end
        if tNow <= prof.tStart
            curCNoDB = settings.defaultCNoDB;
        elseif tNow >= prof.tEnd
            switch shape
                case 'halfsine_dip'
                    curCNoDB = settings.defaultCNoDB;
                otherwise % 'linear_hold'
                    curCNoDB = prof.lowDB;
            end
        else
            frac = (tNow - prof.tStart) / (prof.tEnd - prof.tStart);
            switch shape
                case 'halfsine_dip'
                    curCNoDB = settings.defaultCNoDB - ...
                        (settings.defaultCNoDB - prof.lowDB) * sin(pi*frac);
                otherwise % 'linear_hold'
                    curCNoDB = settings.defaultCNoDB + ...
                        (prof.lowDB - settings.defaultCNoDB) * frac;
            end
        end
        curCNoValuesdB = CNoSetting(satList, curCNoDB) + 1.5;
        curCNoValues   = 10.^(curCNoValuesdB/10);
        carrAmp = 2 * settings.gwnAmp * sqrt(curCNoValues/settings.samplingFreq);
        if settings.fileType == 2
            carrAmp = carrAmp/sqrt(2);
        end
        % Once-per-second progress print, so the profile's actual effect
        % is easy to verify while a run is in progress (and afterward, by
        % scrolling the command window) instead of only showing up as a
        % silent difference in the finished IF file.
        if mod(loopCnt, 1000) == 1
            fprintf('  t=%6.2fs: C/No profile (%s) -> %.2f dB-Hz\n', tNow, shape, curCNoDB);
        end
    end

    for svIndex = 1:length(satList)
        % Compute the transmitting time coresponding to the Rx time -------
        PRN = satList(svIndex);
        [TxTime,satClkErr,satPosBlk] = GetTravelTime(RxTime,RxPosEcef,eph(PRN),settings);
        % include tgd, clock error and relativistic effect (atmospheric
        % delay is applied further below, once elevation/azimuth are known)
        TxTime = TxTime + satClkErr;
        if isNavIC
            % content time stamped in IRNWT, and the receiver's L5 channel
            % delay (a later arrival = an earlier content stamp)
            TxTime = TxTime + navicDtIG(TxTime) - navicRxL5Delay;
        end

        % Elevation/azimuth of this satellite as seen from the receiver's
        % local horizon - needed for antenna gain (if enabled) and always
        % needed for the ionospheric/tropospheric delay below, so it is
        % computed unconditionally. Held constant over this 1 ms block -
        % both antenna gain and atmospheric delay vary far slower than
        % the code/carrier, so per-sample resolution is unnecessary.
        [azRx, elRx] = topocent(RxPosEcef(1,:)', (satPosBlk(1,:) - RxPosEcef(1,:))');

        antGainDb = 0;
        if settings.rxAntennaEnable
            [azRxBody, elRxBody] = applyAttitudeToBodyFrame(azRx, elRx, rollDeg, pitchDeg, headingDeg);
            if isfield(settings,'rxBelowHorizonSlopeDb')
                antGainDb = antGainDb + lookupAntennaGain(rxAntPattern, azRxBody, elRxBody, ...
                    settings.rxBelowHorizonSlopeDb, settings.rxBelowHorizonFloorDb);
            else
                antGainDb = antGainDb + lookupAntennaGain(rxAntPattern, azRxBody, elRxBody);
            end
        end
        if settings.txAntennaEnable
            toEarthCenter = -satPosBlk(1,:);
            toReceiver = RxPosEcef(1,:) - satPosBlk(1,:);
            nadirAngle = acosd(dot(toEarthCenter,toReceiver) / (norm(toEarthCenter)*norm(toReceiver)));
            antGainDb = antGainDb + lookupAntennaGain(txAntPattern, 0, nadirAngle);
        end
        antGainLinear = 10^(antGainDb/20);

        % Ionospheric (Klobuchar) and tropospheric (Saastamoinen) delay.
        % Both are REAL extra propagation delay (the atmosphere always
        % slows the signal down), so - unlike satClkErr, which is a
        % clock-content-stamp offset, not a travel-time change - they
        % push the effective content timestamp EARLIER (subtracted
        % below). Applied uniformly to code and carrier together; the
        % ionosphere's true dispersive code-delay/carrier-advance split
        % is not modeled here.
        atmDelaySec = 0;
        if settings.ionoEnable
            atmDelaySec = atmDelaySec + ionoScale*klobuchar(ionoAlpha, ionoBeta, RxTime(1), ...
                trajectory(loopCnt,1), trajectory(loopCnt,2), elRx, azRx);
        end
        if settings.tropoEnable
            atmDelaySec = atmDelaySec + saastamoinenTropo(elRx, trajectory(loopCnt,1), ...
                trajectory(loopCnt,3), settings.tropoP0, settings.tropoT0, settings.tropoRH);
        end
        TxTime = TxTime - atmDelaySec;

        % generate local code, carrier and Nav data -----------------------
        sapcing = (TxTime(2) - TxTime(1))/blockSize;
        TxTimeSample = linspace(TxTime(1),TxTime(2) - sapcing,blockSize);
        sapcing = (RxTime(2) - RxTime(1))/blockSize;
        RxTimeSample = linspace(RxTime(1),RxTime(2) - sapcing,blockSize);
        % local carrier samples
        if  settings.fileType == 1
            localCarr = cos(carrFreqRad * TxTimeSample - localOsFreq * RxTimeSample);
        elseif settings.fileType == 2
            localCarr = exp(-1i*(carrFreqRad * TxTimeSample - localOsFreq * RxTimeSample));
        end
        
        % Local code samples
        % mod() (not rem()) so the index stays in [1, codeLength] even
        % when TxTimeSample is negative (happens whenever the ephemeris
        % epoch falls exactly at a GPS week boundary, toc==0).
        codePhase = TxTimeSample * settings.codeFreqBasis;
        codeIndex = ceil(mod(codePhase,settings.codeLength)) + 1;
        localCode = caCodeTable(svIndex,codeIndex);
        % Nav bit samples
        if isNavIC
            % NavIC symbols (50 sps): regenerate the 48 s master frame
            % when the transmit time runs past it (each sample is looked
            % up in the frame that contains it, so a block that straddles
            % a frame boundary gets the symbols of both)
            while TxTime(end) >= navBits(svIndex).refTime + 48
                if navBits(svIndex).refTime == -Inf
                    seedTime = TxTime(1);
                else
                    seedTime = navBits(svIndex).refTime + 48 + 0.5;
                end
                [sym, newRef] = navic_generateNavMsg(navBits(svIndex).template, seedTime);
                navBits(svIndex).prevBit = navBits(svIndex).dataBit;
                navBits(svIndex).prevRef = navBits(svIndex).refTime;
                navBits(svIndex).dataBit = reshape(sym.', 1, 2400);
                navBits(svIndex).refTime = newRef;
            end
            symIdx = floor((TxTimeSample - navBits(svIndex).refTime)/bitPeriod) + 1;
            localDataBit = zeros(1, blockSize);
            inCur = symIdx >= 1;
            localDataBit(inCur) = navBits(svIndex).dataBit(symIdx(inCur));
            if any(~inCur)
                prevIdx = floor((TxTimeSample(~inCur) - navBits(svIndex).prevRef)/bitPeriod) + 1;
                localDataBit(~inCur) = navBits(svIndex).prevBit(prevIdx);
            end
        else
        dataTime = TxTimeSample - (navBits(svIndex).refTime - 6);
        bitPhase = ceil(dataTime/bitPeriod);
        localDataBit = navBits(svIndex).dataBit(bitPhase);
        end

        % Update Nav date bits --------------------------------------------
        if ~isNavIC && (TxTime(end) + 1) > (navBits(svIndex).refTime + 30)
            % Time to generate Nav Messages
            navBits(svIndex).ephTime.second = navBits(svIndex).refTime + 30.5;
                        
            % Modulate WN, TOW and CRC into Nav messages
            [navBits(svIndex).dataWord, navBits(svIndex).refTime] = ...
                generateNavMsg(navBits(svIndex).frameMsg,navBits(svIndex).ephTime, ...
                0,navBits(svIndex).dataWord);
            
            % Extract Nav data bits from nav message
            for bitindex = 1:60
                navBits(svIndex).dataBit((bitindex-1) * 30 + 1 : bitindex * 30) = ...
                double(bitget(navBits(svIndex).dataWord(bitindex),30:-1:1)) * 2 - 1;
            end
        end
        
        % Combine each componests of local signals
        localSig = (localCode .* localDataBit .* localCarr) * carrAmp(svIndex) * antGainLinear;
        localSigSum = localSigSum + localSig;

        % Add multipath rays -----------------------------------------------
        % Each ray is the LOS code/nav-bit replica delayed by the ray's
        % extra propagation delay, attenuated, and carrier-rotated by an
        % independently-evolving fading phase Phi_k(t) (Eq. 1 of the
        % referenced paper). The nav bit is reused from the LOS path since
        % the multipath delays here (sub-microsecond to a few chips) are
        % far shorter than the 20 ms nav bit period.
        if settings.multipathEnable
            for k = 1:settings.numMultipathRays
                TxTimeSampleMP = TxTimeSample - mpDelaySec(svIndex,k);

                codePhaseMP = TxTimeSampleMP * settings.codeFreqBasis;
                codeIndexMP = ceil(mod(codePhaseMP,settings.codeLength)) + 1;
                localCodeMP = caCodeTable(svIndex,codeIndexMP);

                fadePhase = 2*pi*mpFadeFreq(svIndex,k)*RxTimeSample + mpPhase0(svIndex,k);
                if settings.fileType == 1
                    localCarrMP = cos(carrFreqRad*TxTimeSampleMP - localOsFreq*RxTimeSample + fadePhase);
                elseif settings.fileType == 2
                    localCarrMP = exp(-1i*(carrFreqRad*TxTimeSampleMP - localOsFreq*RxTimeSample) + 1i*fadePhase);
                end

                localSigMP = (localCodeMP .* localDataBit .* localCarrMP) ...
                    * carrAmp(svIndex) * mpAtt(svIndex,k) * antGainLinear;
                localSigSum = localSigSum + localSigMP;
            end
        end
    end % svIndex = 1:length(satList)

    % Add chirp (sweep) jamming interference -------------------------------
    % Satellite-independent (jams the whole receiver, not one PRN), so
    % added once per block rather than once per satellite. tLocal is
    % measured against absolute simulation time so the sweep phase stays
    % continuous within each sweep period and resets cleanly at tj
    % boundaries, matching a real jammer's VCO retrace.
    if settings.chirpEnable
        blockSpacing = (RxTime(2) - RxTime(1))/blockSize;
        blockTimeSamples = linspace(RxTime(1), RxTime(2) - blockSpacing, blockSize);
        % Elapsed simulation time (seconds since startTime), for gating
        % the jammer to its configured injection window.
        elapsedTimeSamples = blockTimeSamples - startTime;
        injMask = elapsedTimeSamples >= settings.chirpInjStartSec & ...
                  elapsedTimeSamples <  settings.chirpInjEndSec;
        if any(injMask)
            tLocal = mod(blockTimeSamples, settings.chirpSweepPeriodMs/1000);
            jamPhase = 2*pi*(settings.chirpDeltaF0*tLocal + 0.5*chirpK*tLocal.^2) + chirpPhi0;
            if settings.fileType == 1
                localJam = jamAmp * cos(jamPhase);
            elseif settings.fileType == 2
                localJam = jamAmp * exp(-1i*jamPhase);
            end
            localJam(~injMask) = 0;
            localSigSum = localSigSum + localJam;
        end
    end

    % Update Rx time and receiver position --------------------------------
    RxTime(1) = RxTime(2);
    RxTime(2) = RxTime(2) + blockTime;
    % Corresponding receiver positions 
    % next block spans trajectory rows loopCnt+1 .. loopCnt+2 (was loopCnt+1
    % for the end point, which made the signal lag the trajectory by 1 ms)
    RxPosEcef(1,:)  = RxPosEcef(2,:);
    RxPosEcef(2,:)  = trajectory(min(loopCnt + 2, iterCnt + 1),4:6);
    
    % Add noise and filter --------------------------------------------------
    % WGN is added unconditionally - it models receiver thermal noise, not
    % the front-end filter, so disabling the filter (settings.filterEn=0)
    % should not silently produce a noiseless IF file.
    if  settings.fileType == 1
        gwnSamlple = wgn(1,blockSize,settings.gwnAmp^2,'linear');
    elseif settings.fileType == 2
        gwnSamlple = wgn(1,blockSize,settings.gwnAmp^2,'linear')...
            + 1i*wgn(1,blockSize,settings.gwnAmp^2,'linear');
    end
    localSigSum = localSigSum + gwnSamlple;
    if settings.filterEn == 1
        % Do FE filtering, carrying the filter state across blocks
        [localSigSum, filterState] = filter(coef,1,localSigSum,filterState);
    end
    
    % ADC quantization ----------------------------------------------------
    if strcmp(settings.dataType,'int8')
        localSigSum = localSigSum/max(real(localSigSum)) * 2^7;
        quantizedSig = int8(localSigSum);
    elseif strcmp(settings.dataType,'int16')
        localSigSum = localSigSum/max(real(localSigSum)) * 2^12;
        quantizedSig = int16(localSigSum);
    end

    % Write samples. For fileType==2 only the interleaved I/Q write below
    % must run - fwrite() on a complex array silently drops the imaginary
    % part, so writing quantizedSig directly here would corrupt the file.
    if settings.fileType == 1
        fwrite(fid,quantizedSig,settings.dataType);
    elseif settings.fileType == 2
        quantizedSig1 = reshape([real(quantizedSig);imag(quantizedSig)],[],1);
        fwrite(fid,quantizedSig1,settings.dataType);
    end
    
end
%% clear environment
fclose(fid);
close(hwb)
end

function [gammaP, scaleP] = mpDelayParams(elevDeg)
% Gamma-distribution parameters (Table 4) for multipath delay [m] by
% satellite elevation bucket. Wang, Chen & Liu, Sensors 2018, 18, 1149.
gammaTable = [2.62 2.77 2.81 2.56 2.47 2.40];
scaleTable = [129.83 105.52 80.93 65.12 53.22 43.24];
idx = min(max(floor(elevDeg/15) + 1, 1), 6);
gammaP = gammaTable(idx);
scaleP = scaleTable(idx);
end

function lambdaP = mpFadeLambda(elevDeg)
% Rate parameter (Table 8, MEO column) of the double-sided exponential
% distribution of multipath fading frequency [Hz] by elevation bucket.
% Wang, Chen & Liu, Sensors 2018, 18, 1149.
lambdaTable = [17.2 12.3 10.9 10.6 11.2 12.8];
idx = min(max(floor(elevDeg/15) + 1, 1), 6);
lambdaP = lambdaTable(idx);
end

