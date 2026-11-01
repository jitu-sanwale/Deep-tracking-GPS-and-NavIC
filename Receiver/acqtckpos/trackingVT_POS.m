function [TckResultVT, navSolutionsVT] = trackingVT_POS(file,signal,track,cmn, Acquired,cnslxyz,eph,sbf,TckResultCT,navSolutionsCT)
%Purpose:
%   Vector tracking and positioning
%Inputs: 
%	file        - parameters related to the data file to be processed,a structure 
%	signal   	- parameters related to signals,a structure
%	track     	- parameters related to signal tracking,a structure
%	cmn         - parameters commmonly used,a structure
%	Acquired 	- acquisition results
%	cnslxyz 	- initial position in ECEF coordinate
%	eph         - ephemeris 
%	sbf         - parameters used for pseudorange estimation
%	TckResultCT             - conventional tracking results
%	navSolutionsCT          - navigation solutions in conventional tracking mode
%Outputs:
%	TckResultVT         - vector tracking results
%	navSolutionsVT   	- navigation solutions in vector tracking mode
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU


%% initialization
Spacing = [-track.CorrelatorSpacing 0 track.CorrelatorSpacing];
datalength = track.msToProcessVT; 
f0  = signal.codeFreqBasis;
fL1 = signal.Fc;
fs  = signal.Fs;

pdi = track.pdi; % unit:ms integration time
t   = 1e-3;
sv  = Acquired.sv;
svlength    = length(Acquired.sv);
sv_clk      = zeros(1,200);   % GPS 1-32, NavIC 101-114
sv_clk_pos  = zeros(1,200);
eph_idx     = ones(1,svlength);

% Kalman Filter Parameter (state 9: GPS-NavIC inter-system bias, m)
num_state = 9;

% error state vector
error_state = zeros(num_state,1);

% system transition matrix
Dynamic_Model = diag(zeros(1,9));
Dynamic_Model(1,4)  = 1;
Dynamic_Model(2,5)  = 1;
Dynamic_Model(3,6)  = 1;
Dynamic_Model(7,8)  = 1;
Transistion_Matrix  = eye(length(error_state)) + Dynamic_Model * pdi * t;

% error covariance matrix
state_cov = blkdiag(1e5*diag([1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e0,1e0]), 0);   % ISB held at the CT value

% process (System) noise noise covariance matrix
process_noise(1:3,1:3) = diag(ones(1,3)*2e-1);
process_noise(4:6,4:6) = diag(ones(1,3)*1e-1);
process_noise(7,7) = 1e-1;
process_noise(8,8) = 1e-2;
process_noise(9,9) = track.isbPSD*pdi*t;

% measurement noise covariance matrix
mesurement_noise(1:svlength,1:svlength) = eye(svlength)*3e-1;
mesurement_noise(svlength+1:2*svlength,svlength+1:2*svlength) = eye(svlength)*1e-1;

% parameters for measurement noise variance update
flag_corrCovEst2 = 1;
counterUptR = 0;
counter_r = 0;
thresUptR = 200/pdi;

% initialize navigation solutions of VT using conventional positioning resutls
estPos      = navSolutionsCT.usrPos(file.skiptimeVT,:); 
estVel      = navSolutionsCT.usrVel(file.skiptimeVT,:);
clkBias     = navSolutionsCT.clkBias(file.skiptimeVT);
clkDrift    = navSolutionsCT.clkDrift(file.skiptimeVT); 

% total state vector
isb         = 0;
if isfield(navSolutionsCT,'isb'), isb = navSolutionsCT.isb(file.skiptimeVT); end
total_state = [estPos,estVel,clkBias,clkDrift,isb]';

% parameters for tracking loops
carrError       = zeros(1,svlength);    % carrier phase discriminator output
oldCarrError    = zeros(1,svlength);
codeError       = zeros(1,svlength);    % code discriminator output
oldCarrNco      = zeros(1,svlength); 

% parameters for C/N0 estimate
snrIndex	= ones(1,svlength);
K           = 20;
flag_snr    = ones(1,svlength); % flag to calculate C/N0
index_int   = zeros(1,svlength);

% parameter for updating the iono and tropo correction.
corrUpdateSec   = 0.1;
corrUpt         = corrUpdateSec/(pdi*t);
counter_corr    = corrUpt-1 * ones(svlength,1);

% initialize tracking variables of VT using conventional tracking resutls
for svindex = 1 : length(Acquired.sv)
   	prn = sv(svindex);    
    ch                      = gnssChan(prn, signal, file);   % GPS L1 / NavIC L5
    fidCh(svindex)          = ch.fid;
    lamCh(svindex)          = ch.lam;
    ionoSc(svindex)         = ch.ionoScale;
    isN(svindex)            = double(ch.isNavIC);
    codetemp                = ch.code;
    Code(svindex,:)         = [codetemp(end) repmat(codetemp,1,pdi) codetemp(1)];
    codeFreq(svindex)       = TckResultCT(Acquired.sv(svindex)).codeFreq(file.skiptimeVT);
    carrFreq(svindex)       = TckResultCT(Acquired.sv(svindex)).carrFreq(file.skiptimeVT);
    oldCarrFreq(svindex)    = TckResultCT(Acquired.sv(svindex)).carrFreq(file.skiptimeVT);
    remChip(svindex)        = TckResultCT(Acquired.sv(svindex)).remChip(file.skiptimeVT);
    codePhaseStep(svindex)  = codeFreq(svindex)/fs;
    remCarrPhase(svindex)   = TckResultCT(Acquired.sv(svindex)).remCarrPhase(file.skiptimeVT);
    file_ptr(svindex)       = TckResultCT(Acquired.sv(svindex)).absoluteSample(file.skiptimeVT);
    carrFreqBasis(svindex)  = TckResultCT(Acquired.sv(svindex)).carrFreq(file.skiptimeVT);
    codedelay_tck(svindex)  = TckResultCT(Acquired.sv(svindex)).codedelay(file.skiptimeVT);
    oldCarrError(svindex)   = TckResultCT(Acquired.sv(svindex)).carrError(file.skiptimeVT);
    oldCarrNco(svindex)     = TckResultCT(Acquired.sv(svindex)).carrFreq(file.skiptimeVT) - carrFreqBasis(svindex);  
end

if all(isN), isN(:) = 0; end   % NavIC-only: no inter-system bias

% Calculate filter coefficient values
[tau1carr, tau2carr] = calcLoopCoef(track.PLLBW,track.PLLDamp,track.PLLGain);       
bytesPerSample = file.dataPrecision*file.dataType;
fprintf('Vector tracking (VDLL) + EKF positioning, %d ms ...\n', datalength);
%% Start processing
for msIndex = 1: datalength/pdi  
    for svindex = 1 : length(Acquired.sv)
        prn = sv(svindex);        
        
        % read raw data file
        numSample(svindex) = round((signal.codelength*pdi-remChip(svindex)) /(codeFreq(svindex)/fs));
        fseek(fidCh(svindex), file_ptr(svindex)*1,'bof');  	
        rawsignal = readIFData(fidCh(svindex), numSample(svindex), file);
        
        time = eph(prn).TOW1(eph_idx(svindex)) - sbf.sfb1(prn)*1/50 - sbf.nav1(prn)*t ...
                    + msIndex*pdi*t + file.skiptimeVT*t;   
        
        % time to subframe 1, in unit of samples, used for tracking
        diff_of_dat_tck(svindex) = sbf.sfb1(prn)*20*fs*t ...
                                        + sbf.nav1(prn)*fs*t ...
                                        + codedelay_tck(svindex); 
        % time of transmition 
        tot_est_tck(svindex) = eph(prn).TOW1(eph_idx(svindex))...
                                    - diff_of_dat_tck(svindex)/fs...
                                    + msIndex*pdi*t ......
                                    + file.skiptimeVT*t ...
                                    - (1/cmn.cSpeed)*(sv_clk(prn)); % t = t_sv - dt_sv (was +, unlike CT)
                                            
        % Current SV position and Velocity for tck
        [svxyz_tck(svindex,:), sv_vel(svindex,:), sv_clk(prn), sv_clk_vel(prn),  grpdel(prn)] = ...
                                svPosVel(prn,eph,tot_est_tck(svindex),eph_idx(svindex));
        
        %% Iono, trop correction 
        counter_corr(svindex) = counter_corr(svindex) + 1;
        if counter_corr(svindex) ==  corrUpt 
            svenu           = xyz2enu(svxyz_tck(svindex,:),estPos(1:3));
            el_rad(svindex) = atan(svenu(3)/norm(svenu(1:2)));
            az_rad(svindex) = (pi/2)-atan2(svenu(1),svenu(2));
            az(svindex)     = az_rad(svindex) * 180/pi;
            el(svindex)     = el_rad(svindex) * 180/pi;
            
            temp = xyz2llh(estPos);
            user_ll	= [temp(1:2) .* 180 /pi temp(3)];
            ionodel(svindex)        = ionoSc(svindex)*ionocorr(tot_est_tck(svindex),svxyz_tck(svindex,:), estPos(1:3));
            tropodel_unb3(svindex)  = abs(trop_UNB3(cmn.doy,user_ll(1),user_ll(3),el(svindex)));
            counter_corr(svindex)   = 0;
        end
        
        
        %% Predict code freq
        r = sqrt(sum((svxyz_tck(svindex,:) - estPos(1:3)).^2));
        predictedPr_tck(svindex) = r + clkBias + isN(svindex)*isb - sv_clk(prn) + grpdel(prn)*cmn.cSpeed + tropodel_unb3(svindex) + ionodel(svindex);  % predicted raw pseudorange, Groves (14.142); signs as in trackingUT_POS prCorr
        % Earth rotation over the geometric travel time (passing the pseudorange,
        % which carries the ~7 ms receiver clock bias, rotated the SV by ~14 m)
        svxyzr_tck(svindex,:) = erotcorr(svxyz_tck(svindex,:),r);
        
        r = sqrt(sum((svxyzr_tck(svindex,:) - estPos(1:3)).^2));
        predictedPr_tck(svindex) = r + clkBias + isN(svindex)*isb - sv_clk(prn) + grpdel(prn)*cmn.cSpeed + tropodel_unb3(svindex) + ionodel(svindex);
        
        if msIndex == 1
            codeFreq(svindex) = TckResultCT(Acquired.sv(svindex)).codeFreq(file.skiptimeVT+1);  
        else
            codeFreq(svindex) = f0*(1-(predictedPr_tck(svindex)-predictedPr_last(svindex))/(cmn.cSpeed*pdi*t));
            % a diverged navigation solution must not drive the code NCO outside
            % any physical code Doppler (|f - f0| < 10 kHz covers > 2.9 km/s)
            if ~isfinite(codeFreq(svindex)) || abs(codeFreq(svindex) - f0) > 1e4
                codeFreq(svindex) = f0 + sign(codeFreq(svindex) - f0)*1e4*isfinite(codeFreq(svindex));
            end             
        end
        predictedPr_last(svindex) = predictedPr_tck(svindex);
        
        codePhaseStep(svindex) = codeFreq(svindex)/fs ;
         
        t_CodeEarly      = (0 + Spacing(1) + remChip(svindex)) : codePhaseStep(svindex) : ...
                                ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(1) + remChip(svindex));
        t_CodePrompt     = (0 + Spacing(2) + remChip(svindex)) : codePhaseStep(svindex) : ...
                                ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(2) + remChip(svindex));
        t_CodeLate       = (0 + Spacing(3) + remChip(svindex)) : codePhaseStep(svindex) : ...
                                ((numSample(svindex) -1) * codePhaseStep(svindex) + Spacing(3) + remChip(svindex));
        
        CodeEarly 	= Code(svindex,(ceil(t_CodeEarly) + 1));
        CodePrompt	= Code(svindex,(ceil(t_CodePrompt) + 1));
        CodeLate	= Code(svindex,(ceil(t_CodeLate) + 1));
        
        remChip(svindex) = (t_CodePrompt(numSample(svindex)) + codePhaseStep(svindex)) - signal.codelength*pdi;  
                            
        CarrTime = (0: numSample(svindex))./signal.Fs;
        Wave = 2*pi*((carrFreq(svindex)).* CarrTime) + remCarrPhase(svindex);
        remCarrPhase(svindex) = rem(Wave(numSample(svindex)+1), 2*pi);
        [InphaseSignal, QuadratureSignal] = carrierWipeoff(rawsignal, Wave(1:numSample(svindex)), signal);
        
        E_i	= sum(CodeEarly    .*InphaseSignal);  
        E_q = sum(CodeEarly    .*QuadratureSignal);
        P_i	= sum(CodePrompt   .*InphaseSignal);  
        P_q = sum(CodePrompt   .*QuadratureSignal);
        L_i	= sum(CodeLate     .*InphaseSignal);  
        L_q = sum(CodeLate     .*QuadratureSignal); 
        
        % Calculate CN0
        if (flag_snr(svindex) == 1)
            index_int(svindex) = index_int(svindex) + 1;
            Zk(svindex,index_int(svindex)) = P_i^2 + P_q^2;
            if mod(index_int(svindex),K) == 0
                CN0_VT(snrIndex(svindex),svindex) = cn0Moment(Zk(svindex,:), 1*t*pdi, track.cn0Floor);   % floor for an invalid estimate
                index_int(svindex)  = 0;
                snrIndex(svindex)   = snrIndex(svindex) + 1;
            end
        end
        
        % PLL discriminator
        carrError(svindex)      = atan(P_q/P_i)/(2.0*pi);
        carrNco(svindex)        = oldCarrNco(svindex) + (tau2carr/tau1carr) * (carrError(svindex) - ...
                                        oldCarrError(svindex)) + carrError(svindex) * (pdi*1e-3/tau1carr);
        oldCarrNco(svindex)     = carrNco(svindex);
        oldCarrError(svindex)   = carrError(svindex);
        carrFreq(svindex)       = carrFreqBasis(svindex) + carrNco(svindex);
        oldCarrFreq(svindex)    = carrFreq(svindex);
        
        % DLL discriminator
        E	= sqrt(E_i^2+E_q^2);
        L	= sqrt(L_i^2+L_q^2); 
        codeError(svindex)  = 0.5*(E-L)/(E+L);
        
        % form the pseudorange error measurement
        Z(msIndex,svindex) = (codeError(svindex))*cmn.cSpeed / codeFreq(svindex);
        delta_Z_ch(msIndex,svindex) = Z(msIndex,svindex); 
        
        % Result Record
        TckResultVT(prn).P_i(msIndex)                = P_i;
        TckResultVT(prn).P_q(msIndex)                = P_q;
        TckResultVT(prn).E_i(msIndex)                = E_i;
        TckResultVT(prn).E_q(msIndex)                = E_q;
        TckResultVT(prn).L_i(msIndex)                = L_i;
        TckResultVT(prn).L_q(msIndex)                = L_q;
        TckResultVT(prn).carrError(msIndex)          = carrError(svindex);
        TckResultVT(prn).codeError(msIndex)          = codeError(svindex);
        TckResultVT(prn).codePhase(msIndex)          = remChip(svindex);
        TckResultVT(prn).carrPhase(msIndex)          = remCarrPhase(svindex);
        TckResultVT(prn).codeFreq(msIndex)           = codeFreq(svindex);
        TckResultVT(prn).carrFreq(msIndex)           = carrFreq(svindex);
        TckResultVT(prn).absoluteSample(msIndex)     = ftell(fidCh(svindex)); 
        file_ptr(svindex)                            = TckResultVT(prn).absoluteSample(msIndex);
        TckResultVT(prn).sv_vel(msIndex,:)           = sv_vel(svindex,:); 
        cdMod = mod(TckResultVT(prn).absoluteSample(msIndex)/bytesPerSample - remChip(svindex)/codePhaseStep(svindex), fs*t); % sub-sample code epoch
        TckResultVT(prn).codedelay(msIndex)          = unwrapCodeDelay(codedelay_tck(svindex), cdMod, fs*t);
        
        codedelay_tck(svindex)                       = TckResultVT(prn).codedelay(msIndex);
        
    end % end for svindex in Tracking
    
    %% interpolate for positioning
    if msIndex == 1
        for svindex = 1 : svlength
            prn = sv(svindex);
            x1 = [TckResultCT(prn).absoluteSample(file.skiptimeVT-1)/(file.dataPrecision*file.dataType), ... 
                        (TckResultVT(prn).absoluteSample(msIndex))/(file.dataPrecision*file.dataType)];
            y1 = [0, delta_Z_ch(msIndex,svindex)];
            y2 = [TckResultCT(Acquired.sv(svindex)).codedelay(round(file.skiptimeVT-1)), TckResultVT(prn).codedelay(msIndex)];
            y3 = [TckResultCT(Acquired.sv(svindex)).carrFreq(round(file.skiptimeVT-1)), carrFreq(svindex)];     
            x_s = (file.skip + file.skiptimeVT)*fs*t + msIndex*pdi*fs*t;
            Z(msIndex,svindex)         = safeInterp(x1,y1,x_s);
            codedelay_pos(svindex)     = safeInterp(x1,y2,x_s);
            new_carrFreq(svindex)      = safeInterp(x1,y3,x_s);
            
            sv_clk_pos = sv_clk;            
            
            oldabsoluteSample_pos(svindex) = (TckResultVT(prn).absoluteSample(msIndex));
        end
    else
        for svindex = 1 : svlength
            prn = sv(svindex);
            x1 = [oldabsoluteSample_pos(svindex)/(file.dataPrecision*file.dataType), ...
                        TckResultVT(prn).absoluteSample(msIndex)/(file.dataPrecision*file.dataType)]; 
            y1 = [oldcode(svindex), delta_Z_ch(msIndex,svindex)] ;
            y2 = [oldcodedelay_pos(svindex), TckResultVT(prn).codedelay(msIndex)];
            y3 = [oldnew_carrFreq(svindex), carrFreq(svindex)];
            x_s = (file.skip + file.skiptimeVT)*fs*t + msIndex*pdi*fs*t;
            Z(msIndex,svindex)         = safeInterp(x1,y1,x_s);
            codedelay_pos(svindex)     = safeInterp(x1,y2,x_s);
            new_carrFreq(svindex)      = safeInterp(x1,y3,x_s);
        
            oldabsoluteSample_pos(svindex) = TckResultVT(prn).absoluteSample(msIndex);
        end
    end
    
    oldcode = Z(msIndex,1:svlength);
    oldcodedelay_pos = codedelay_pos;
    oldnew_carrFreq =  new_carrFreq;
     
    for svindex = 1 : svlength
        prn = sv(svindex);
        
        diff_of_dat_pos(svindex) = sbf.sfb1(prn)*20*fs*t ...
                                        + sbf.nav1(prn)*fs*t ...
                                        + (codedelay_pos(svindex));                

        tot_est_pos(svindex) = eph(prn).TOW1(eph_idx(svindex))...
                                    - diff_of_dat_pos(svindex)/fs ...
                                    + msIndex*pdi*t ...
                                    + file.skiptimeVT*t ...
                                    - (1/cmn.cSpeed)*sv_clk_pos(prn); % t = t_sv - dt_sv (was +, unlike CT)
     
        [svxyz_pos(svindex,:), sv_vel_pos(svindex,:), sv_clk_pos(prn), sv_clk_vel(prn), grpdel(prn)] = ...
                                    svPosVel(prn,eph,tot_est_pos(svindex),eph_idx(svindex));
        
        r = sqrt(sum((svxyz_pos(svindex,:) - estPos(1:3)).^2));
        predictedPr_pos(svindex) = r  + clkBias + isN(svindex)*isb - sv_clk_pos(prn) + grpdel(prn)*cmn.cSpeed ...
                                        + tropodel_unb3(svindex) + ionodel(svindex);    % 
        svxyzr_pos(svindex,:) = erotcorr(svxyz_pos(svindex,:),r); % geometric travel time, see above
        r = sqrt(sum((svxyzr_pos(svindex,:) - estPos(1:3)).^2));
        a_pos(svindex,:) = (svxyzr_pos(svindex,:)-estPos(1:3)) / r;
        % the GPS-NavIC inter-system bias is carried over from CT (not
        % re-estimated): the VDLL code-discriminator measurements carry
        % per-channel offsets of a few metres, whose NavIC mean the bias state
        % would otherwise absorb at the handover transient
        H_pos(svindex,:) = [-a_pos(svindex,:) 0 0 0 1 0 0];
        H_pos(svindex+svlength,:) = [0 0 0 -a_pos(svindex,:) 0 1 0];

        % measured pseudorange rate
        prr_measured(svindex)	= (new_carrFreq(svindex) - signal.IF)*lamCh(svindex);
        
        % predicted pseudorange rate 
        prr_predicted(svindex)	= (estVel - sv_vel_pos(svindex,:))*a_pos(svindex,:)';

        % pseudorange rate error measurement
        Z(msIndex,svlength+svindex) = (prr_predicted(svindex) - prr_measured(svindex)) - clkDrift + sv_clk_vel(prn);    
    end
    
    newZ = Z(msIndex,:); % measurements
    
    error_state = zeros(num_state,1);
    error_state = Transistion_Matrix * error_state;       
    
    state_cov = Transistion_Matrix * state_cov * transpose(Transistion_Matrix) + process_noise; % error covariance of the predicted state 
    % channel screening (C/N0, elevation); VT keeps tracking excluded channels
    % from the navigation solution, so they can be re-acquired when the signal returns
    cn0Now = 99*ones(1,svlength);
    for svindex = 1:svlength
        if snrIndex(svindex) > 1, cn0Now(svindex) = CN0_VT(snrIndex(svindex)-1,svindex); end
    end
    [rows, usedCh] = measMask(cn0Now, el, [], msIndex, track);
    navSolutionsVT.usedCh(msIndex,:) = usedCh;
    kalman_gain = zeros(num_state, 2*svlength);
    kalman_gain(:,rows) = state_cov * transpose(H_pos(rows,:)) * inv(H_pos(rows,:) * state_cov * transpose(H_pos(rows,:)) + mesurement_noise(rows,rows)); % Kalman gain   
        
    counterUptR = counterUptR + 1;  % counter for update measurement noise variance, R
    recordR(counterUptR,:) = ((newZ' - H_pos * error_state));    
    
    error_state = error_state + kalman_gain * (newZ' - H_pos * error_state); % error state update   
    state_cov = (eye(num_state) - kalman_gain * H_pos) * state_cov; % error covariance update
        
    total_state = total_state + error_state; % total state update
    estPos = total_state(1:3)';
    estVel = total_state(4:6)';
    clkBias = total_state(7);
    clkDrift = total_state(8);
    isb      = total_state(9);
    
    %% record navigation results
    llh     = xyz2llh(estPos); 
    L_b     = llh(1);
    lamda_b = llh(2); 
    C_e_n = [ -sin(lamda_b)           cos(lamda_b)         	 0;...
              -sin(L_b)*cos(lamda_b) -sin(L_b)*sin(lamda_b)	 cos(L_b);...
              -cos(L_b)*cos(lamda_b) -cos(L_b)*sin(lamda_b)	-sin(L_b);];    
    usrenuvel(msIndex,:) = C_e_n * estVel';
    
    usrenu(msIndex,:) = xyz2enu(estPos,cnslxyz);  
    usrllh(msIndex,:) = xyz2llh(estPos);
    usrllh(msIndex,1:2)	= usrllh(msIndex,1:2)*180/pi;   
    
    navSolutionsVT.usrTime(msIndex,:)        = time;  
    navSolutionsVT.usrPos(msIndex,:)         = estPos;
    navSolutionsVT.usrVel(msIndex,:)         = estVel;
    navSolutionsVT.usrPosENU(msIndex,:)      = usrenu(msIndex,:);
    navSolutionsVT.usrVelENU(msIndex,:)      = usrenuvel(msIndex,:);
    navSolutionsVT.usrPosLLH(msIndex,:)   	 = usrllh(msIndex,:);
    navSolutionsVT.clkDrift(msIndex,:)  	 = [time,clkDrift];
    navSolutionsVT.clkBias(msIndex,:)        = [time,clkBias];
    navSolutionsVT.isb(msIndex,1)            = isb;
    navSolutionsVT.state(msIndex,:)          = error_state;
    navSolutionsVT.svxyz_pos(:,:,msIndex)    = svxyzr_pos;
    navSolutionsVT.svvel_pos(:,:,msIndex)    = sv_vel_pos;
    navSolutionsVT.kalman_gain(:,:,msIndex)	 = kalman_gain;
    navSolutionsVT.state_cov(msIndex,:)      = diag(state_cov);
    navSolutionsVT.meas_inno(msIndex,:)      = ((newZ' - H_pos * error_state));
    navSolutionsVT.newZ(msIndex,:)           = newZ;
    navSolutionsVT.predicted_z(msIndex,:)    = H_pos * error_state;
    navSolutionsVT.satEA(msIndex,:)      = el;
    navSolutionsVT.satAZ(msIndex,:)      = az;
    
    
    %% predict postion and clkBias at next epoch
    total_state = Transistion_Matrix * total_state;
    estPos = total_state(1:3)';
    clkBias = total_state(7)'; 
    
    
    %% update R by measurement variance 
    if flag_corrCovEst2 == 1 && counterUptR == thresUptR
        tmpR = diag(var(recordR));
        mesurement_noise(1:svlength,1:svlength) = ...
                                                    tmpR(1:svlength,1:svlength);     
        mesurement_noise(svlength+1:2*svlength,svlength+1:2*svlength) = ...
                                                    tmpR(svlength+1:2*svlength,svlength+1:2*svlength); 
        for idx = 1 : svlength
            if mesurement_noise(idx,idx) >= 12000
                mesurement_noise(idx,idx) = 12000;
            elseif mesurement_noise(idx,idx) <= 0.01	
                mesurement_noise(idx,idx) = 0.01;	
            end
            if mesurement_noise(idx+svlength,idx+svlength) >= 400
                mesurement_noise(idx+svlength,idx+svlength) = 400;
            elseif mesurement_noise(idx+svlength,idx+svlength) <= 0.01	
                mesurement_noise(idx+svlength,idx+svlength) = 0.01;	
            end
        end
        counterUptR = 0;
        counter_r = counter_r + 1;
        navSolutionsVT.R(counter_r,:) = diag(mesurement_noise);
    end   
    
    %%
    TOW_USR(msIndex) = time;
    if mod(msIndex, track.progressEvery) == 0
        fprintf('VT: index = %6d  %6.4f E %+9.3f N %+9.3f U %+9.3f  B %f  D %+3.4f\n',...
            msIndex,time,usrenu(msIndex,1),usrenu(msIndex,2),usrenu(msIndex,3),clkBias, clkDrift);
    end
   
end % end for msIndex 

save(['navSolVT_',file.fileName,'.mat'], 'navSolutionsVT','eph','TOW_USR');
save(['tckRstVT_',file.fileName,'.mat'], 'TckResultVT','CN0_VT');


