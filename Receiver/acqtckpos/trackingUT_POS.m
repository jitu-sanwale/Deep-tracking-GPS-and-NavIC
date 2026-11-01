function [TckResultUT, navSolutionsUT] = trackingUT_POS(file,signal,track,cmn,Acquired,cnslxyz,eph,sbf,TckResultCT,navSolutionsCT,imu,aiding)
%Purpose:
%   Ultratight GPS/INS tracking and positioning - VDCaP architecture
%   (Sun, Petovello & Cannon, "Ultratight GPS/Reduced-IMU Integration for
%   Land Vehicle Navigation", IEEE TAES 49(3), 2013, Sec. II-B), here with a
%   full 6-DOF IMU:
%     code    - VDLL: the code NCO follows the pseudorange rate predicted from
%               the INS and is steered onto the pseudorange predicted by the
%               navigation filter; the discriminator is a measurement only.
%     carrier - CaP (cascaded / Doppler-aided PLL, Fig. 3): carrier NCO =
%               navigation-derived Doppler f_Nav + 2nd-order loop filter
%               output f_LF (Eq. 4). Doppler measurement f_dm = f_Nav + f_LF
%               (Eq. 3, f_LF averaged over the update interval). After every
%               filter update the loop filter is reset when PLI < 0.9.
%     filter  - 18-state error-state EKF, closed loop on an ECEF strapdown
%               INS (insNavECEF.m): attitude(3) velocity(3) position(3)
%               accel bias(3) gyro bias(3) clock bias(1) clock drift(1)
%               (Groves, 2nd ed., Sec. 14.2 tightly coupled ECEF model),
%               plus the GPS-NavIC inter-system bias (state 18, m), seen by
%               the NavIC L5 pseudoranges only (18 states with NavIC).
%               The INS runs at the IMU rate, the filter is updated every
%               track.utUpdateMs ms with pseudoranges and pseudorange rates.
%   Started from the CT solution at file.skiptimeVT, like trackingVT_POS.
%   aiding = 'none' gives the vector-tracking receiver (VT) on the same
%   framework: identical VDLL and CaP channels, measurements, noise models,
%   gating and feedback, with the INS replaced by a 12-state Singer
%   position-velocity-acceleration model [vel pos acc clk drift isb]
%   (Singer, IEEE TAES 6(4), 1970). The only difference is the INS aiding.
%Inputs:
%	file, signal, track, cmn, Acquired, eph, sbf  - as trackingVT_POS
%	cnslxyz         - ENU reference position (ECEF)
%	TckResultCT     - conventional tracking results (handover)
%	navSolutionsCT  - conventional navigation solution (handover, prBase0)
%	imu             - IMU struct from generateIMU_iNAT_FSLG01.m: t [s from
%	                  the start of the data file], f, w (body-frame interval
%	                  averages), truth.attRad (for the initial attitude only);
%	                  [] for aiding 'none'
%	aiding          - 'INS' (ultratight, default) or 'none' (vector tracking)
%Outputs:
%	TckResultUT     - tracking results
%	navSolutionsUT  - navigation solutions at the filter update epochs

%% initialization
if nargin < 12 || isempty(aiding), aiding = 'INS'; end
ins = strcmpi(aiding, 'INS');
if ins, tag = 'UT'; else, tag = 'VT'; end
capTransfer = isfield(track, 'capNavTransfer') && track.capNavTransfer.(tag);
Spacing  = [-track.CorrelatorSpacing 0 track.CorrelatorSpacing];
datalength = track.msToProcessVT;
f0  = signal.codeFreqBasis;
fL1 = signal.Fc;
fs  = signal.Fs;
c   = cmn.cSpeed;
lamL1   = c/fL1;
lamCode = c/f0;
T   = 1e-3;                       % coherent integration (pdi = 1 ms)
L   = track.utUpdateMs;           % filter update interval [ms]
sv  = Acquired.sv;
svlength = length(sv);
sv_clk   = zeros(1,200);  sv_clk_vel = zeros(1,200);  grpdel = zeros(1,200);   % GPS 1-32, NavIC 101-114
eph_idx  = ones(1,svlength);
bytesPerSample = file.dataPrecision*file.dataType;
omgIE = 7.292115e-5;
OmgIE = [0 -omgIE 0; omgIE 0 0; 0 0 0];
bSec  = 75e-3;                    % as pr_est.m
base0 = navSolutionsCT.prBase0;
tOff  = (file.skip + file.skiptimeVT)*T;   % data-file time of msIndex = 0

%% INS / filter initialisation from the CT solution at the handover epoch
k0 = file.skiptimeVT;
rIns  = navSolutionsCT.usrPos(k0,:)';
vIns  = navSolutionsCT.usrVel(k0,:)';
clkB  = navSolutionsCT.clkBias(k0);
clkD  = navSolutionsCT.clkDrift(k0);
isb   = 0;
if isfield(navSolutionsCT,'isb'), isb = navSolutionsCT.isb(k0); end
tIns  = tOff;                     % data-file time of the INS state
ba = zeros(3,1);  bg = zeros(3,1);
aK = zeros(3,1);                  % VT: kinematic ECEF acceleration estimate
d2r = pi/180; g0 = 9.80665;
Sclk = track.utClkPhasePSD;  Sdrf = track.utClkFreqPSD;
Sisb = track.isbPSD;
aLast = zeros(3,1);
if ins
% Initial attitude: handed over at the iNAT-FSLG-01 datasheet GNSS-aided
% accuracy (roll/pitch 0.01 deg, heading 0.02 deg rms) - applied as a fixed
% error on the truth attitude, since a few seconds of flight are not enough
% for an in-flight alignment and attitude is not part of the CT solution.
attT  = interp1(imu.truth.t, imu.truth.attRad, tIns, 'linear', 'extrap');
att0  = attT + track.utAttErr0Deg*pi/180;
llh0  = xyz2llh(rIns');
Cbe   = cen(llh0(1), llh0(2))' * cnbEuler(att0)';

P = diag([ (track.utP0AttDeg*d2r)^2*ones(1,3), track.utP0Vel^2*ones(1,3), track.utP0Pos^2*ones(1,3), ...
           (300e-6*g0)^2*ones(1,3), (0.05*d2r/3600)^2*ones(1,3), track.utP0Clk^2, track.utP0Drift^2, 1 ]);
nx = 18;                                     % 17 INS/clock states + inter-system bias
iV = 4:6;  iR = 7:9;  iB = 16;  iD = 17;  iI = 18;   % state indices shared with VT

% process noise PSDs (datasheet; unmodelled scale factor / misalignment /
% non-linearity of the accelerometers added as white noise, see qAccel below)
ds = imu.datasheet;
Sg   = ds.gyroARW^2;                         % rad^2/s
Sa   = ds.accVRW^2;                          % m^2/s^3
Sba  = 2*ds.accBI^2/ds.tauBI;                % Gauss-Markov driving PSD
Sbg  = 2*ds.gyroBI^2/ds.tauBI;

kImu  = find(imu.t > tIns + 1e-9, 1);        % next IMU sample to integrate
else
% VT: error states [vel pos acc] (estimated minus true) + clock bias, drift
% and ISB (true minus estimated), as the UT states they replace
Cbe = eye(3);
P  = diag([ track.utP0Vel^2*ones(1,3), track.utP0Pos^2*ones(1,3), track.vtP0Acc^2*ones(1,3), ...
            track.utP0Clk^2, track.utP0Drift^2, 1 ]);
nx = 12;
iV = 1:3;  iR = 4:6;  iA = 7:9;  iB = 10;  iD = 11;  iI = 12;
tauA = track.vtSingerTau;
Sacc = 2*track.vtSingerSigmaA^2/tauA;        % Singer driving PSD [m^2/s^5]
end

%% tracking-loop initialisation (handover from CT, as trackingVT_POS)
if ins, capBW = track.utPLLBW; else, capBW = track.vtPLLBW; end   % residual CaP bandwidth [Hz]
[wn, zeta] = deal(capBW*8*0.707/(4*0.707^2+1), 0.707);
fLF   = zeros(1,svlength);  fSum = zeros(1,svlength);
codeSteer = zeros(1,svlength);
for svindex = 1:svlength
    prn = sv(svindex);
    ch                     = gnssChan(prn, signal, file);   % GPS L1 / NavIC L5
    fidCh(svindex)         = ch.fid;
    lamCh(svindex)         = ch.lam;
    ionoSc(svindex)        = ch.ionoScale;
    isN(svindex)           = double(ch.isNavIC);
    codetemp               = ch.code;
    Code(svindex,:)        = [codetemp(end) codetemp codetemp(1)];
    codeFreq(svindex)      = TckResultCT(prn).codeFreq(k0);
    carrFreq(svindex)      = TckResultCT(prn).carrFreq(k0);
    remChip(svindex)       = TckResultCT(prn).remChip(k0);
    remCarrPhase(svindex)  = TckResultCT(prn).remCarrPhase(k0);
    file_ptr(svindex)      = TckResultCT(prn).absoluteSample(k0);
    codedelay_tck(svindex) = TckResultCT(prn).codedelay(k0);
end
if all(isN), isN(:) = 0; end   % NavIC-only: no inter-system bias
el = navSolutionsCT.satEA(k0,:);  az = navSolutionsCT.satAZ(k0,:);
ionodel = zeros(1,svlength); tropodel = zeros(1,svlength);
corrUpt = 100;  counter_corr = (corrUpt-1)*ones(1,svlength);
firstNco = true(1,svlength);    % fSum is seeded from CT's carrier at the first ms

% C/N0 estimation (as trackingVT_POS)
snrIndex = ones(1,svlength); K = 20; index_int = zeros(1,svlength);
Zk = zeros(svlength,K); CN0_UT = zeros(1,svlength);

% interval accumulators
accCode = zeros(1,svlength); accFLF = zeros(1,svlength);
accNBD  = zeros(1,svlength); accNBP = zeros(1,svlength); accN = 0;
cohI    = zeros(1,svlength); cohQ   = zeros(1,svlength);

nEpoch = floor(datalength/L);
navSolutionsUT.usrTime   = zeros(nEpoch,1);
navSolutionsUT.tData     = zeros(nEpoch,1);
navSolutionsUT.usrPos    = zeros(nEpoch,3);  navSolutionsUT.usrVel    = zeros(nEpoch,3);
navSolutionsUT.usrPosENU = zeros(nEpoch,3);  navSolutionsUT.usrVelENU = zeros(nEpoch,3);
navSolutionsUT.usrPosLLH = zeros(nEpoch,3);  navSolutionsUT.attEuler  = zeros(nEpoch,3);
navSolutionsUT.ba = zeros(nEpoch,3); navSolutionsUT.bg = zeros(nEpoch,3);
navSolutionsUT.clkBias = zeros(nEpoch,1); navSolutionsUT.clkDrift = zeros(nEpoch,1);
navSolutionsUT.stateStd  = zeros(nEpoch,nx);
navSolutionsUT.isb       = zeros(nEpoch,1);
navSolutionsUT.usedCh    = false(nEpoch,svlength);
navSolutionsUT.PLI       = zeros(nEpoch,svlength);
navSolutionsUT.reset     = false(nEpoch,svlength);
navSolutionsUT.innovPr   = zeros(nEpoch,svlength); navSolutionsUT.innovPrr = zeros(nEpoch,svlength);
navSolutionsUT.CN0       = zeros(nEpoch,svlength);

if ins
    fprintf('Ultratight VDCaP (VDLL + CaP, 18-state EKF + 6-DOF INS), %d ms, update every %d ms ...\n', datalength, L);
else
    fprintf('Vector tracking VDCaP (VDLL + CaP, 12-state Singer PVA EKF, no INS), %d ms, update every %d ms ...\n', datalength, L);
end
nE = 0;
%% Start processing
for msIndex = 1:datalength
    tEpoch = tOff + msIndex*T;          % data-file time of this epoch's sample x_s

    % INS: integrate every IMU sample up to the epoch (bias-corrected) and
    % propagate the error covariance
    if ~ins
        % VT: Singer PVA propagation to the epoch (no inertial measurements)
        tau = tEpoch - tIns;
        if tau > 0
            rIns = rIns + vIns*tau + 0.5*aK*tau^2;
            vIns = vIns + aK*tau;
            aK   = exp(-tau/tauA)*aK;
            clkB = clkB + clkD*tau;
            [Phi, Q] = pvaSystemModel(tau, tauA, Sacc, Sclk, Sdrf, Sisb);
            P = Phi*P*Phi' + Q;
            tIns = tEpoch;  aLast = aK;
        end
    end
    while ins && kImu <= numel(imu.t) && imu.t(kImu) <= tEpoch + 1e-9
        tau = imu.t(kImu) - tIns;
        fb  = imu.f(kImu,:)' - ba;   wb = imu.w(kImu,:)' - bg;
        CbeOld = Cbe;  rOld = rIns;
        [rIns, vIns, Cbe] = insNavECEF(tau, rIns, vIns, Cbe, fb, wb);
        clkB = clkB + clkD*tau;
        fie  = CbeOld*fb;
        aLast = fie + gravityECEF(rIns) - 2*OmgIE*vIns;
        [Phi, Q] = utSystemModel(tau, CbeOld, fie, rOld, OmgIE, Sg, Sa + qAccel(fb, ds), Sba, Sbg, Sclk, Sdrf, Sisb);
        P = Phi*P*Phi' + Q;
        tIns = imu.t(kImu);  kImu = kImu + 1;
    end

    for svindex = 1:svlength
        prn = sv(svindex);

        % read raw data
        numSample(svindex) = round((signal.codelength - remChip(svindex)) / (codeFreq(svindex)/fs));
        fseek(fidCh(svindex), file_ptr(svindex), 'bof');
        rawsignal = readIFData(fidCh(svindex), numSample(svindex), file);

        % satellite state at transmission (as trackingVT_POS, t = t_sv - dt_sv)
        diff_of_dat_tck = sbf.sfb1(prn)*20*fs*T + sbf.nav1(prn)*fs*T + codedelay_tck(svindex);
        tot_tck = eph(prn).TOW1(eph_idx(svindex)) - diff_of_dat_tck/fs + msIndex*T + file.skiptimeVT*T ...
                    - sv_clk(prn)/c;
        [svp, svv, sv_clk(prn), sv_clk_vel(prn), grpdel(prn)] = svPosVel(prn,eph,tot_tck,eph_idx(svindex));

        % receiver state at the middle of this correlation block (INS extrapolated)
        tMid = (file_ptr(svindex)/bytesPerSample + numSample(svindex)/2)/fs;
        dtx  = tMid - tIns;
        rU   = rIns + vIns*dtx;   vU = vIns + aLast*dtx;   dU = clkD;
        rgeo = norm(svp' - rU);
        svpr = erotcorr(svp, rgeo)';   % erotcorr returns a column
        u    = (svpr' - rU)/norm(svpr' - rU);
        prrNav = u'*(svv' - vU) + dU - sv_clk_vel(prn);     % pseudorange rate incl. clocks [m/s]
        fNav   = -prrNav/lamCh(svindex);                     % navigation-derived Doppler [Hz] (L1 or L5)

        % iono / tropo for the predicted pseudorange (every 100 ms)
        counter_corr(svindex) = counter_corr(svindex) + 1;
        if counter_corr(svindex) == corrUpt
            svenu = xyz2enu(svpr, rU');
            el(svindex) = atan(svenu(3)/norm(svenu(1:2)))*180/pi;
            az(svindex) = ((pi/2)-atan2(svenu(1),svenu(2)))*180/pi;
            ull = xyz2llh(rU');
            ionodel(svindex)  = ionoSc(svindex)*ionocorr(tot_tck, svpr, rU');   % L1 Klobuchar -> L5 for NavIC
            tropodel(svindex) = rxTropo(ull(1)*180/pi, ull(3), el(svindex));   % same model as the simulator
            counter_corr(svindex) = 0;
        end

        % NCOs: VDLL code rate from the INS, steered onto the predicted pseudorange;
        % CaP carrier = navigation Doppler + loop filter
        codeFreq(svindex) = f0*(1 - prrNav/c) + codeSteer(svindex);
        if firstNco(svindex)
            fSum(svindex) = carrFreq(svindex) - signal.IF - fNav;   % continue CT's carrier
            fLF(svindex)  = fSum(svindex);
            firstNco(svindex) = false;
        end
        carrFreq(svindex) = signal.IF + fNav + fLF(svindex);
        codePhaseStep = codeFreq(svindex)/fs;

        tP = (0:numSample(svindex)-1)*codePhaseStep + remChip(svindex);
        CodeEarly  = Code(svindex, ceil(tP + Spacing(1)) + 1);
        CodePrompt = Code(svindex, ceil(tP + Spacing(2)) + 1);
        CodeLate   = Code(svindex, ceil(tP + Spacing(3)) + 1);
        remChip(svindex) = tP(end) + codePhaseStep - signal.codelength;

        CarrTime = (0:numSample(svindex))./fs;
        Wave = 2*pi*(carrFreq(svindex).*CarrTime) + remCarrPhase(svindex);
        remCarrPhase(svindex) = rem(Wave(numSample(svindex)+1), 2*pi);
        [I, Qs] = carrierWipeoff(rawsignal, Wave(1:numSample(svindex)), signal);

        E_i = sum(CodeEarly.*I);  E_q = sum(CodeEarly.*Qs);
        P_i = sum(CodePrompt.*I); P_q = sum(CodePrompt.*Qs);
        L_i = sum(CodeLate.*I);   L_q = sum(CodeLate.*Qs);

        % C/N0
        index_int(svindex) = index_int(svindex) + 1;
        Zk(svindex,index_int(svindex)) = P_i^2 + P_q^2;
        if index_int(svindex) == K
            CN0_UT(snrIndex(svindex),svindex) = cn0Moment(Zk(svindex,:), T, track.cn0Floor);   % floor for an invalid estimate
            index_int(svindex) = 0; snrIndex(svindex) = snrIndex(svindex) + 1;
        end

        % CaP: Costas discriminator [cycles] and 2nd-order loop filter (Eq. 4)
        carrError = atan(P_q/P_i)/(2*pi);
        fSum(svindex) = fSum(svindex) + wn^2*T*carrError;
        fLF(svindex)  = fSum(svindex) + 2*zeta*wn*carrError;

        % VDLL discriminator (measurement only)
        E = sqrt(E_i^2+E_q^2);  Lm = sqrt(L_i^2+L_q^2);
        codeError = 0.5*(E-Lm)/(E+Lm);

        accCode(svindex) = accCode(svindex) + codeError;
        accFLF(svindex)  = accFLF(svindex) + fLF(svindex);
        % PLI from utPLIcohMs-ms coherent sums: with 1-ms powers the noise
        % biases (I^2-Q^2)/(I^2+Q^2) to ~0.85 at 38 dB-Hz even when locked
        cohI(svindex) = cohI(svindex) + P_i;  cohQ(svindex) = cohQ(svindex) + P_q;
        if mod(msIndex, track.utPLIcohMs) == 0
            accNBD(svindex) = accNBD(svindex) + cohI(svindex)^2 - cohQ(svindex)^2;
            accNBP(svindex) = accNBP(svindex) + cohI(svindex)^2 + cohQ(svindex)^2;
            cohI(svindex) = 0;  cohQ(svindex) = 0;
        end

        % record
        TckResultUT(prn).P_i(msIndex) = P_i;  TckResultUT(prn).P_q(msIndex) = P_q;
        TckResultUT(prn).E_i(msIndex) = E_i;  TckResultUT(prn).E_q(msIndex) = E_q;
        TckResultUT(prn).L_i(msIndex) = L_i;  TckResultUT(prn).L_q(msIndex) = L_q;
        TckResultUT(prn).carrError(msIndex) = carrError;
        TckResultUT(prn).codeError(msIndex) = codeError;
        TckResultUT(prn).codeFreq(msIndex)  = codeFreq(svindex);
        TckResultUT(prn).carrFreq(msIndex)  = carrFreq(svindex);
        TckResultUT(prn).fNav(msIndex)      = fNav;
        TckResultUT(prn).fLF(msIndex)       = fLF(svindex);
        TckResultUT(prn).absoluteSample(msIndex) = ftell(fidCh(svindex));
        file_ptr(svindex) = TckResultUT(prn).absoluteSample(msIndex);
        cdMod = mod(file_ptr(svindex)/bytesPerSample - remChip(svindex)/codePhaseStep, fs*T);
        codedelay_new = unwrapCodeDelay(codedelay_tck(svindex), cdMod, fs*T);
        TckResultUT(prn).codedelay(msIndex) = codedelay_new;

        % code delay interpolated to the common epoch sample x_s
        xs = (file.skip + file.skiptimeVT + msIndex)*fs*T;
        if msIndex == 1
            xPrev = TckResultCT(prn).absoluteSample(k0)/bytesPerSample;
        else
            xPrev = TckResultUT(prn).absoluteSample(msIndex-1)/bytesPerSample;
        end
        codedelay_pos(svindex) = safeInterp([xPrev, file_ptr(svindex)/bytesPerSample], ...
                                            [codedelay_tck(svindex), codedelay_new], xs);
        codedelay_tck(svindex) = codedelay_new;
    end
    accN = accN + 1;

    %% filter update every L ms
    if mod(msIndex, L) ~= 0, continue; end
    nE = nE + 1;
    dtx = tEpoch - tIns;                          % 0 when the IMU grid hits the epoch
    rU = rIns + vIns*dtx;  vU = vIns + aLast*dtx;  bU = clkB + clkD*dtx;

    cn0Now = 99*ones(1,svlength);
    for svindex = 1:svlength
        if snrIndex(svindex) > 1, cn0Now(svindex) = CN0_UT(snrIndex(svindex)-1,svindex); end
    end

    H = zeros(2*svlength,nx);  dz = zeros(2*svlength,1);  R = zeros(2*svlength);
    uAll = zeros(3,svlength);
    prPred = zeros(1,svlength); prRep = zeros(1,svlength); svState = zeros(svlength,7); prCorr = zeros(1,svlength);
    for svindex = 1:svlength
        prn = sv(svindex);
        diff_of_dat_pos = sbf.sfb1(prn)*20*fs*T + sbf.nav1(prn)*fs*T + codedelay_pos(svindex);
        tot_pos = eph(prn).TOW1(eph_idx(svindex)) - diff_of_dat_pos/fs + msIndex*T + file.skiptimeVT*T ...
                    - sv_clk(prn)/c;
        [svp, svv, sv_clk(prn), sv_clk_vel(prn), grpdel(prn)] = svPosVel(prn,eph,tot_pos,eph_idx(svindex));
        svpr = erotcorr(svp, norm(svp' - rU))';
        rng_ = norm(svpr' - rU);  u = (svpr' - rU)/rng_;
        svState(svindex,:) = [svpr, svv, sv_clk_vel(prn)];
        prCorr(svindex)    = - sv_clk(prn) + grpdel(prn)*c + ionodel(svindex) + tropodel(svindex);

        % replica pseudorange (receiver time origin of pr_est / CT) and the
        % VDLL measurement = replica + averaged discriminator
        prRep(svindex)  = c*(bSec + (sbf.sfb1(prn)*20*fs*T + sbf.nav1(prn)*fs*T + codedelay_pos(svindex) - 1 - base0)/fs);
        prMeas          = prRep(svindex) + (accCode(svindex)/accN)*lamCode;
        prPred(svindex) = rng_ + bU + isN(svindex)*isb + prCorr(svindex);

        % CaP Doppler measurement: f_dm = f_Nav + <f_LF>  ->  prr innovation = -lambda*<f_LF>
        prrPred = u'*(svv' - vU) + clkD - sv_clk_vel(prn);
        prrMeas = prrPred - lamCh(svindex)*(accFLF(svindex)/accN);

        dz(svindex)          = prMeas - prPred(svindex);
        dz(svindex+svlength) = prrMeas - prrPred;
        uAll(:,svindex) = u;
        H(svindex, iR) = u';             H(svindex, iB) = 1;   H(svindex, iI) = isN(svindex);
        H(svindex+svlength, iV) = u';    H(svindex+svlength, iD) = 1;

        cn0 = 10^(cn0Now(svindex)/10);
        R(svindex,svindex) = lamCode^2*(1/(4*T*cn0))*(1 + 2/(T*cn0))/L + track.utPrFloor^2;   % d = 1 chip
        R(svindex+svlength,svindex+svlength) = track.utPrrSigma^2*(10^4.5/cn0)*(1 + 1/(2*T*cn0)) + 0.02^2;
    end

    % channel selection: C/N0 and elevation mask (measMask) + innovation gate
    [rows, usedCh] = measMask(cn0Now, el, [], msIndex, track);
    S  = H*P*H' + R;
    nis = abs(dz)./sqrt(diag(S));
    if msIndex > track.gateStart
        bad = (nis(1:svlength) > track.utGate) | (nis(svlength+1:end) > track.utGate);
        keep = usedCh & ~bad';
        if sum(keep) >= 4, usedCh = keep; end
        rows = [find(usedCh), find(usedCh)+svlength];
    end

    Hs = H(rows,:);  Rs = R(rows,rows);
    Kg = P*Hs'/(Hs*P*Hs' + Rs);
    x  = Kg*dz(rows);
    IKH = eye(nx) - Kg*Hs;
    P  = IKH*P*IKH' + Kg*Rs*Kg';

    % closed-loop correction of the INS, IMU biases and clock
    if ins
        Cbe  = (eye(3) - skew(x(1:3)))*Cbe;
        ba = ba + x(10:12);    bg = bg + x(13:15);
    else
        aK = aK - x(iA);  aLast = aK;
    end
    vIns = vIns - x(iV);  rIns = rIns - x(iR);
    clkB = clkB + x(iB);   clkD = clkD + x(iD);   isb = isb + x(iI);
    if capTransfer
        % the update changes the predicted range rate by u'*dv + d(drift);
        % hand that Doppler over from the loop filter to f_Nav
        dfNav = -(uAll'*x(iV) + x(iD))'./lamCh;
        fSum = fSum - dfNav;  fLF = fLF - dfNav;
    end
    rU = rIns + vIns*dtx;  vU = vIns + aLast*dtx;  bU = clkB + clkD*dtx;

    % CaP filter reset (PLI < 0.9 after the update) and VDLL steering
    PLI = accNBD./max(accNBP, eps);
    resetCh = PLI < track.utPLIreset;
    fSum(resetCh) = 0;  fLF(resetCh) = 0;
    for svindex = 1:svlength
        % steer the replica onto the corrected predicted pseudorange over utCodeSteerTau
        prPredNew = norm(svState(svindex,1:3)' - rU) + bU + isN(svindex)*isb + prCorr(svindex);
        codeSteer(svindex) = (prRep(svindex) - prPredNew)/(lamCode*track.utCodeSteerTau);
    end

    % record
    llh  = xyz2llh(rU');
    Cen  = cen(llh(1), llh(2));
    navSolutionsUT.usrTime(nE)      = eph(sv(1)).TOW1(1) - sbf.sfb1(sv(1))/50 - sbf.nav1(sv(1))*T + (msIndex + file.skiptimeVT)*T;
    navSolutionsUT.tData(nE)        = tEpoch;
    navSolutionsUT.msIndex(nE,1)    = msIndex;
    navSolutionsUT.usrPos(nE,:)     = rU';
    navSolutionsUT.usrVel(nE,:)     = vU';
    navSolutionsUT.usrPosENU(nE,:)  = xyz2enu(rU', cnslxyz);
    navSolutionsUT.usrVelENU(nE,:)  = ([0 1 0; 1 0 0; 0 0 -1]*Cen*vU)';
    navSolutionsUT.usrPosLLH(nE,:)  = [llh(1:2)*180/pi, llh(3)];
    if ins, navSolutionsUT.attEuler(nE,:) = ctmToEuler((Cen*Cbe)')*180/pi; end   % roll pitch yaw (deg)
    navSolutionsUT.ba(nE,:) = ba';  navSolutionsUT.bg(nE,:) = bg';
    navSolutionsUT.clkBias(nE)      = bU;
    navSolutionsUT.clkDrift(nE)     = clkD;
    navSolutionsUT.isb(nE)          = isb;
    navSolutionsUT.stateStd(nE,:)   = sqrt(diag(P))';
    navSolutionsUT.usedCh(nE,:)     = usedCh;
    navSolutionsUT.PLI(nE,:)        = PLI;
    navSolutionsUT.reset(nE,:)      = resetCh;
    navSolutionsUT.innovPr(nE,:)    = dz(1:svlength)';
    navSolutionsUT.innovPrr(nE,:)   = dz(svlength+1:end)';
    navSolutionsUT.CN0(nE,:)        = cn0Now;
    navSolutionsUT.satEA(nE,:)      = el;
    navSolutionsUT.satAZ(nE,:)      = az;

    accCode(:) = 0; accFLF(:) = 0; accNBD(:) = 0; accNBP(:) = 0; accN = 0;

    if mod(msIndex, track.progressEvery) == 0
        enu = navSolutionsUT.usrPosENU(nE,:);
        fprintf('%s: index = %6d  E %+9.3f N %+9.3f U %+9.3f  B %.3f  D %+.4f  used %d/%d  resets %d\n', ...
            tag, msIndex, enu(1), enu(2), enu(3), bU, clkD, sum(usedCh), svlength, sum(resetCh));
    end
end
navSolutionsUT.refPos = cnslxyz;

if ins
    save(['navSolUT_',file.fileName,'.mat'], 'navSolutionsUT','eph');
    save(['tckRstUT_',file.fileName,'.mat'], 'TckResultUT','CN0_UT');
else
    navSolutionsVT = navSolutionsUT;  TckResultVT = TckResultUT;  CN0_VT = CN0_UT; %#ok<NASGU>
    save(['navSolVT_',file.fileName,'.mat'], 'navSolutionsVT','eph');
    save(['tckRstVT_',file.fileName,'.mat'], 'TckResultVT','CN0_VT');
end
end

%% ------------------------------------------------------------------------
function [Phi, Q] = pvaSystemModel(tau, tauA, Sacc, Sclk, Sdrf, Sisb)
% VT: Singer position-velocity-acceleration model in ECEF, error states
% [vel pos acc clk drift isb]; first-order transition as utSystemModel.
Phi = eye(12);
Phi(1:3,7:9) = eye(3)*tau;                 % velocity <- acceleration
Phi(4:6,1:3) = eye(3)*tau;                 % position <- velocity
Phi(7:9,7:9) = eye(3)*(1 - tau/tauA);      % acceleration: first-order Gauss-Markov
Phi(10,11)   = tau;
Q = diag([zeros(1,6), Sacc*ones(1,3), Sclk, Sdrf, Sisb])*tau;
end

%% ------------------------------------------------------------------------
function [Phi, Q] = utSystemModel(tau, Cbe, fie, r, OmgIE, Sg, Sa, Sba, Sbg, Sclk, Sdrf, Sisb)
% Groves (2nd ed.) Eq. 14.48-14.50 / TC_KF_Epoch: ECEF error-state transition
% and process noise, states [att vel pos ba bg clk drift].
R0 = 6378137; e = 0.0818191908425;
llh = xyz2llh(r');  sL = sin(llh(1));
reS = R0/sqrt(1 - e^2*sL^2)*sqrt(cos(llh(1))^2 + (1 - e^2)^2*sL^2);
g0  = norm(gravityECEF(r));
Phi = eye(18);
Phi(1:3,1:3)   = eye(3) - OmgIE*tau;
Phi(1:3,13:15) = Cbe*tau;
Phi(4:6,1:3)   = -tau*skew(fie);
Phi(4:6,4:6)   = eye(3) - 2*OmgIE*tau;
Phi(4:6,7:9)   = tau*2*g0/reS*(r*r')/(r'*r);   % +2g0/reS: Groves (14.49) has the downward vector gamma, i.e. -g0*r/|r|
Phi(4:6,10:12) = Cbe*tau;
Phi(7:9,4:6)   = eye(3)*tau;
Phi(16,17)     = tau;
Q = diag([Sg*ones(1,3), Sa(:)', zeros(1,3), Sba*ones(1,3), Sbg*ones(1,3), Sclk, Sdrf, Sisb])*tau;   % ISB: random walk
end

function S = qAccel(fb, ds)
% Unmodelled accelerometer errors (scale factor, misalignment, quadratic
% non-linearity) treated as white noise with a 1 s correlation time; at the
% ~60 g peak the 50 ug/g^2 non-linearity alone is up to 1.8 m/s^2.
a = abs(fb(:));
S = ((ds.accSF*a).^2 + (ds.accMis*norm(fb))^2 + (ds.accNL*a.^2).^2)*1;
end

function C = cen(lat, lon)
% ECEF -> NED
sL = sin(lat); cL = cos(lat); sl = sin(lon); cl = cos(lon);
C = [-sL*cl -sL*sl cL; -sl cl 0; -cL*cl -cL*sl -sL];
end

function C = cnbEuler(eul)
% Euler (roll pitch yaw) -> C_n^b (Groves Euler_to_CTM)
sp = sin(eul(1)); cp = cos(eul(1)); st = sin(eul(2)); ct = cos(eul(2)); ss = sin(eul(3)); cs = cos(eul(3));
C = [ct*cs, ct*ss, -st; -cp*ss+sp*st*cs, cp*cs+sp*st*ss, sp*ct; sp*ss+cp*st*cs, -sp*cs+cp*st*ss, cp*ct];
end

function eul = ctmToEuler(Cnb)
% C_n^b -> Euler (roll pitch yaw) (Groves CTM_to_Euler)
eul = [atan2(Cnb(2,3), Cnb(3,3)), -asin(max(min(Cnb(1,3),1),-1)), atan2(Cnb(1,2), Cnb(1,1))];
end

function S = skew(a)
S = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
end
