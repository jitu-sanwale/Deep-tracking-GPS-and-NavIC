function [TckResultCT, navSolutionsCT] = trackingCT_POS(file,signal,track,cmn, Acquired,cnslxyz,eph,sbf)
%Purpose:
%   Conventional tracking and positioning using EKF and WLS
%Inputs: 
%	file        - parameters related to the data file to be processed,a structure 
%	signal   	- parameters related to signals,a structure
%	track     	- parameters related to signal tracking,a structure
%	cmn         - parameters commmonly used,a structure
%	Acquired 	- acquisition results
%	cnslxyz 	- initial position in ECEF coordinate
%	eph         - ephemeris
%	sbf         - parameters used for pseudorange estimation
%
%Outputs:
%	TckResultCT         - conventional tracking results
%	navSolutionsCT   	- navigation solutions in conventional tracking mode
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU


%% 
% If no reference position is given, the receiver position and clock are
% initialised with a least-squares fix at the first epoch, and that fix is
% used as the ENU reference.
useLSref = isempty(cnslxyz);
if useLSref, cnslxyz = [0 0 0]; end
sv_clk              = zeros(1,200);   % GPS 1-32, NavIC 101-114
clkBias_kf         	= 0; 
clkBias_wls         	= 0; 
clkDrift            = 0;  
estusr_kf           = cnslxyz;%
estusr_wls           = cnslxyz;%
estVel              = zeros(1,3); 

Spacing = [-track.CorrelatorSpacing 0 track.CorrelatorSpacing];

sv      = Acquired.sv;
f0      = signal.codeFreqBasis;
fL1     = signal.Fc;
fs      = signal.Fs ;
pdi     = track.pdi; 
t       = signal.ms;
svlength    = length(sv);
datalength  = track.msToProcessCT;

% Kalman Filter Parameter: position, velocity, clock bias, clock drift and,
% for GPS + NavIC, the GPS-NavIC inter-system bias (state 9, metres) that
% only the NavIC pseudoranges see
num_state   = 9;

% error state vector
error_state = zeros(num_state,1);
total_state = [cnslxyz,zeros(1,6)]';%zeros(num_state,1);

% system transition matrix
Dynamic_Model = diag(zeros(1,num_state));   % ISB: random walk
Dynamic_Model(1,4)  = 1;
Dynamic_Model(2,5)  = 1;
Dynamic_Model(3,6)  = 1;
Dynamic_Model(7,8)  = 1;
Transistion_Matrix  = eye(length(error_state)) + Dynamic_Model*pdi*t;

% error covariance matrix
state_cov = blkdiag(1e5*diag([1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e0,1e0]), track.isbP0^2);

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

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(track.DLLBW,track.DLLDamp,track.DLLGain);
[tau1carr, tau2carr] = calcLoopCoef(track.PLLBW,track.PLLDamp,track.PLLGain);

% initialize tracking parameters using acquisition results
for svindex = 1:length(sv)
    prn                     = sv(svindex);
    ch                      = gnssChan(prn, signal, file);   % GPS L1 / NavIC L5
    fidCh(svindex)          = ch.fid;
    lamCh(svindex)          = ch.lam;
    ionoSc(svindex)         = ch.ionoScale;
    isN(svindex)            = double(ch.isNavIC);
    codetemp                = ch.code;
    Code(svindex,:)         = [codetemp(end) repmat(codetemp,1,pdi) codetemp(1)];
    AcqCodeDelay(svindex)   = Acquired.codedelay(svindex);
    file_ptr(svindex)       = signal.Sample - AcqCodeDelay(svindex) -1 + file.skip *fs*t;  
    carrFreq(svindex)       = Acquired.fineFreq(svindex);
    AcqFreq(svindex)        = Acquired.fineFreq(svindex);
    
    oldcodedelay_pos(svindex) = 0;
    oldabsoluteSample_pos(svindex) = 0;
end

% NavIC-only receiver: no GPS rows, so the clock is referred to GPST through
% message type 9 alone and the inter-system bias is not estimated
if all(isN), isN(:) = 0; end

% parameters for C/N0 estimate
snrIndex	= ones(1,svlength);
K           = 20;
flag_snr    = ones(1,svlength); % flag to calculate C/N0
index_int   = zeros(1,svlength);

eph_idx     = ones(1,svlength);

corrUpdateSec   = 0.1;
corrUpt         = corrUpdateSec/(pdi*t);
counter_corr    = corrUpt-1 * ones(svlength,1); 
        
% Tracking parameters
carrNco      = zeros(1,svlength);
oldCarrNco  = zeros(1,svlength);
oldCarrError       = zeros(1,svlength);
codeNco         = zeros(1,svlength);
code_outputLast     = zeros(1,svlength);
DLLdiscriLast       = zeros(1,svlength);
remChip             = zeros(1,svlength); 
codeFreq            = ones(1,svlength)*f0; 
remCarrPhase        = zeros(1,svlength);
carrError           = zeros(1,svlength);
codeError           = zeros(1,svlength);
delayValue          = zeros(svlength,datalength/pdi); 
codedelayUnw        = zeros(1,svlength);   % continuous (unwrapped) code delay, samples
bytesPerSample      = file.dataPrecision*file.dataType;
 
fprintf('Conventional tracking + EKF positioning, %d ms ...\n', datalength);
%%
for msIndex = 1: datalength/pdi 
    for svindex = 1 :svlength
        prn = sv(svindex);
        
        % read raw data file
        numSample = round((signal.codelength*pdi-remChip(svindex))/(codeFreq(svindex)/signal.Fs));  
               
        delayValue(svindex,msIndex) = numSample - signal.Sample*pdi;
        
        fseek(fidCh(svindex), file_ptr(svindex)*file.dataPrecision*file.dataType,'bof');  
        rawsignal = readIFData(fidCh(svindex), numSample, file);

        file_ptr(svindex)   = file_ptr(svindex) + numSample;  
        
        t_CodeEarly 	= (0 + Spacing(1) + remChip(svindex)) : codeFreq(svindex)/signal.Fs : ((numSample -1) * ...
                                codeFreq(svindex)/signal.Fs + Spacing(1) + remChip(svindex));
        t_CodePrompt 	= (0 + Spacing(2) + remChip(svindex)) : codeFreq(svindex)/signal.Fs : ((numSample -1) * ...
                                codeFreq(svindex)/signal.Fs + Spacing(2) + remChip(svindex));
        t_CodeLate   	= (0 + Spacing(3) + remChip(svindex)) : codeFreq(svindex)/signal.Fs : ((numSample -1) * ...
                                codeFreq(svindex)/signal.Fs + Spacing(3) + remChip(svindex));
                            
        CodeEarly 	= Code(svindex,(ceil(t_CodeEarly) + 1));
        CodePrompt	= Code(svindex,(ceil(t_CodePrompt) + 1));
        CodeLate	= Code(svindex,(ceil(t_CodeLate) + 1));
        remChip(svindex) = t_CodePrompt(numSample) + codeFreq(svindex)/signal.Fs - signal.codelength*pdi;
        stepUsed = codeFreq(svindex)/signal.Fs;  % code phase step used for this block
        
        CarrTime = (0:numSample)./signal.Fs;
        Wave = 2*pi*((carrFreq(svindex)).*CarrTime) + remCarrPhase(svindex);
        remCarrPhase(svindex) = rem(Wave(numSample+1), 2*pi);
        [InphaseSignal, QuadratureSignal] = carrierWipeoff(rawsignal, Wave(1:numSample), signal);
        
        E_i = sum(CodeEarly    .*InphaseSignal);  
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
                CN0_CT(snrIndex(svindex),svindex) = cn0Moment(Zk(svindex,:), 1*t*pdi, track.cn0Floor);   % floor for an invalid estimate
                index_int(svindex)  = 0;
                snrIndex(svindex)   = snrIndex(svindex) + 1;
            end
        end
        
        % Implement code loop filter and generate NCO command
        E = sqrt(E_i^2+E_q^2);
        L = sqrt(L_i^2+L_q^2);
        codeError(svindex) = 0.5*(E-L)/(E+L);  % DLL discriminator
        codeNco(svindex) = code_outputLast(svindex) + (tau2code/tau1code)*(codeError(svindex)...
                                    - DLLdiscriLast(svindex)) + codeError(svindex)* ((pdi*t)/tau1code);
        DLLdiscriLast(svindex) = codeError(svindex);
        code_outputLast(svindex) = codeNco(svindex);
        codeFreq(svindex) = signal.codeFreqBasis - codeNco(svindex);
        
        % PLL discriminator
        carrError(svindex) = atan(P_q/P_i)/(2*pi);  % PLL discriminator
        carrNco(svindex) = oldCarrNco(svindex) + (tau2carr/tau1carr)*(carrError(svindex) ...
                                        - oldCarrError(svindex)) + carrError(svindex) * ((pdi*t)/tau1carr);
        oldCarrNco(svindex) = carrNco(svindex); 
        oldCarrError(svindex) = carrError(svindex);
        carrFreq(svindex)  = AcqFreq(svindex) + carrNco(svindex);  % Modify carrier freq 
        
        % Data Recording
        TckResultCT(prn).P_i(msIndex)             = P_i;
        TckResultCT(prn).P_q(msIndex)             = P_q;
        TckResultCT(prn).E_i(msIndex)             = E_i;
        TckResultCT(prn).E_q(msIndex)             = E_q;
        TckResultCT(prn).L_i(msIndex)             = L_i;
        TckResultCT(prn).L_q(msIndex)             = L_q;
        TckResultCT(prn).carrError(msIndex)       = carrError(svindex);
        TckResultCT(prn).codeError(msIndex)       = codeError(svindex);
        TckResultCT(prn).codeFreq(msIndex)        = codeFreq(svindex);
        TckResultCT(prn).carrFreq(msIndex)        = carrFreq(svindex);
        TckResultCT(prn).numSample(msIndex)       = numSample;
        TckResultCT(prn).remChip(msIndex)         = remChip(svindex);        
        TckResultCT(prn).remCarrPhase(msIndex)    = remCarrPhase(svindex);
        TckResultCT(prn).absoluteSample(msIndex)  = ftell(fidCh(svindex)); 
        TckResultCT(prn).absoluteSampleCodedelay(msIndex)  = mod(TckResultCT(prn).absoluteSample(msIndex)/...
                                                                (file.dataPrecision*file.dataType),fs*t);
        % code-epoch position with sub-sample resolution: the epoch lies remChip
        % chips before the next sample (1 sample = 60 m at 5 MHz, so integer
        % samples alone give tens of metres of pseudorange quantisation)
        cdMod = mod(TckResultCT(prn).absoluteSample(msIndex)/bytesPerSample - remChip(svindex)/stepUsed, fs*t);
        if msIndex == 1
            codedelayUnw(svindex) = cdMod;
        else
            codedelayUnw(svindex) = unwrapCodeDelay(codedelayUnw(svindex), cdMod, fs*t);
        end
        TckResultCT(prn).codedelay(msIndex)       = codedelayUnw(svindex);
        TckResultCT(prn).delayValue(msIndex)      = delayValue(svindex,msIndex);
        
    end % end for svindex in Tracking
    
	%% Interpolate
    for svindex=1:svlength
    	prn = sv(svindex);
        x1 = [oldabsoluteSample_pos(svindex)/(file.dataPrecision*file.dataType) ...
                TckResultCT(prn).absoluteSample(msIndex)/(file.dataPrecision*file.dataType)];
      	y2 = [oldcodedelay_pos(svindex) TckResultCT(prn).codedelay(msIndex)];
      	x_s = file.skip*fs*t + msIndex*pdi*fs*t ;%+ (sbf.nav1(prn))*fs*t 
        codedelay_pos(msIndex, svindex) = safeInterp(x1,y2,x_s);
        oldabsoluteSample_pos(svindex) = TckResultCT(prn).absoluteSample(msIndex);
    end
    oldcodedelay_pos = codedelay_pos(msIndex, :);
    
   
    %% Pseudorange calculation and correction     
 [pseudorange, ~, prBase0] = pr_est(Acquired, sbf, signal,TckResultCT, msIndex);
 
    %% Least-squares initialisation of position and clock at the first epoch
    if msIndex == 1
        pr0 = zeros(1,svlength); svp0 = zeros(svlength,3);
        for svindex = 1 : svlength
            prn = sv(svindex);
            diff0 = sbf.sfb1(prn)*20*fs*t + sbf.nav1(prn)*fs*t + codedelay_pos(msIndex,svindex);
            tot0  = eph(prn).TOW1(eph_idx(svindex)) - diff0/fs + pdi*msIndex*t;
            [svp0(svindex,:), ~, clk0, ~, gd0] = svPosVel(prn,eph,tot0,eph_idx(svindex));
            pr0(svindex) = pseudorange(svindex) + clk0 - gd0*cmn.cSpeed;
        end
        [x0, b0, okLS, isb0] = lsInitFix(pr0, svp0, isN);
        llh0 = xyz2llh(x0);
        fprintf('LS initial fix: lat %.6f lon %.6f h %.1f m, clock %.1f m (converged=%d)\n', ...
                llh0(1)*180/pi, llh0(2)*180/pi, llh0(3), b0, okLS);
        estusr_kf = x0; clkBias_kf = b0;
        total_state(1:3) = x0'; total_state(7) = b0; total_state(9) = isb0;
        isb_kf = isb0;
        if any(isN), fprintf('LS initial GPS-NavIC inter-system bias %.2f m\n', isb0); end
        if useLSref, cnslxyz = x0; end
    end
 
    usr_clk = clkBias_kf;
    estusr = estusr_kf;
    isb_kf  = total_state(9);
    
    for svindex = 1 : svlength
        prn = sv(svindex);
        
        %
        time = eph(prn).TOW1(eph_idx(svindex)) - sbf.sfb1(prn) * 1/50 ...
              	- sbf.nav1(prn) * 1/1000 + msIndex*pdi*1/1000;
            
        % difference for all satellites in unit of samples        
        fs_new = (f0/codeFreq(svindex))*fs;
        diff_of_dat_pos(svindex) = sbf.sfb1(prn)*20*fs*t ...     
                                        + sbf.nav1(prn)*fs*t ...
                                        + (codedelay_pos(msIndex,svindex));                                    
        % time of transmition                          
        tot_est_pos(svindex) = eph(prn).TOW1(eph_idx(svindex)) ...   % TOW1: TOW of current subframe 1
                                    - diff_of_dat_pos(svindex)/fs ...% 
                                    + pdi*msIndex*t ...
                                    - (1/cmn.cSpeed)*(sv_clk(prn)); 

        [svxyz(svindex,:), sv_vel(svindex,:), sv_clk(prn), sv_clk_vel(prn), grpdel] = ...
                                        svPosVel(prn,eph,tot_est_pos(svindex),eph_idx(svindex));
        
        % C/A-code pseudorange corrected for satellite clock and Tgd 
        prvec(svindex)      = pseudorange(svindex) - usr_clk - isN(svindex)*isb_kf + sv_clk(prn) - grpdel*cmn.cSpeed ;    
        
        % Adjust satellite position coordinates for earth rotation correction
        svxyzr(svindex,:)   = erotcorr(svxyz(svindex,:),prvec(svindex));   
        svenu               = xyz2enu(svxyzr(svindex,:),estusr(1:3));
        
        % tropospheric and ionospheric delay correction	
        counter_corr(svindex) = counter_corr(svindex) + 1;
        if counter_corr(svindex) ==  corrUpt
            svenu           = xyz2enu(svxyzr(svindex,:),estusr(1:3));%
            el_rad(svindex) = atan(svenu(3)/norm(svenu(1:2)));
            az_rad(svindex) = (pi/2)-atan2(svenu(1),svenu(2));
            az(svindex)     = az_rad(svindex)*180/pi;
            el(svindex)     = el_rad(svindex)*180/pi;
            temp            = xyz2llh(estusr(1:3));
            user_ll         = [temp(1:2).*180/pi temp(3)];
            ionodel(svindex)        = ionoSc(svindex)*ionocorr(tot_est_pos(svindex),svxyzr(svindex,:), estusr(1:3));  % L1 Klobuchar -> L5 for NavIC
            tropodel_unb3(svindex)  = abs(trop_UNB3(cmn.doy,user_ll(1),user_ll(3),el(svindex))); 
            counter_corr(svindex)   = 0;     
        end
        
        prvec(svindex) = prvec(svindex) - ionodel(svindex) - tropodel_unb3(svindex); 
    end 

    
    %% Position cal using LS method

    % Position cal using KF
    for svindex = 1 : svlength
        r                       = sqrt(sum((svxyzr(svindex,:) - estusr_kf(1:3)).^2));
        a(svindex,:)            = (svxyzr(svindex,:)-estusr_kf(1:3))/r;
        H(svindex,:)            = [-a(svindex,:) 0 0 0 1 0 isN(svindex)];
        H(svindex+svlength,:)   = [0 0 0 -a(svindex,:) 0 1 0];
        
        d_p(svindex,1)          = (TckResultCT(Acquired.sv(svindex)).carrFreq(msIndex) ...
                                        - signal.IF)*lamCh(svindex) + sum( sv_vel(svindex,:) .* a(svindex,:));
 
        % measured pseudorange rate
        prr_measured(svindex)	= (TckResultCT(Acquired.sv(svindex)).carrFreq(msIndex) - signal.IF)*lamCh(svindex) ;
        
        %predicted pseudorange rate 
        prr_predicted(svindex)	= (estVel-sv_vel(svindex,:))*a(svindex,:)';
       
        d_p(svindex,1) = prr_predicted(svindex) - prr_measured(svindex) - clkDrift + sv_clk_vel(prn);% 
        pr_delta(svindex,1) = prvec(svindex) - norm(svxyzr(svindex,:) - estusr_kf(1:3));% ; 
    end
    
    Z = [pr_delta; d_p]; % measurements
    	
	error_state = zeros(num_state,1);
    error_state = Transistion_Matrix * error_state; % predict
    state_cov = Transistion_Matrix * state_cov * transpose(Transistion_Matrix) + process_noise;
    
    % channel screening (C/N0, elevation, innovation gate)
    cn0Now = 99*ones(1,svlength);
    for svindex = 1:svlength
        if snrIndex(svindex) > 1, cn0Now(svindex) = CN0_CT(snrIndex(svindex)-1,svindex); end
    end
    [rows, usedCh] = measMask(cn0Now, el, pr_delta, msIndex, track);
    navSolutionsKF.usedCh(msIndex,:) = usedCh;
    kalman_gain = zeros(num_state, 2*svlength);
    kalman_gain(:,rows) = state_cov * transpose(H(rows,:)) * inv(H(rows,:) * state_cov * transpose(H(rows,:)) + mesurement_noise(rows,rows));
    
    counterUptR = counterUptR + 1; % counter for update measurement noise variance, R
    recordR(counterUptR,:) = (Z - H*error_state);% measurement innovation
    
    error_state = error_state + kalman_gain*(Z-H*error_state);
    state_cov = (eye(num_state) - kalman_gain * H) * state_cov ;
 
    total_state  = total_state + error_state;
    estusr_kf   = total_state(1:3)';
    estVel      = total_state(4:6)';
    clkBias_kf  = total_state(7);
    clkDrift    = total_state(8);
       
    %record results
 	llh     = xyz2llh(estusr_kf);  
    L_b     = llh(1);
    lamda_b = llh(2);
    C_e_n = [ -sin(lamda_b)           cos(lamda_b)         	 0;...
              -sin(L_b)*cos(lamda_b) -sin(L_b)*sin(lamda_b)	 cos(L_b);...
              -cos(L_b)*cos(lamda_b) -cos(L_b)*sin(lamda_b)	-sin(L_b);];  
    usrenuvel(msIndex,:) =  C_e_n * estVel';
    
    usrenu(msIndex,:)                   = xyz2enu(estusr_kf(1:3),cnslxyz);
    usrllh(msIndex,:)                   = xyz2llh(estusr_kf(1:3));
    usrllh(msIndex,1:2)                 = usrllh(msIndex,1:2)*180/pi;
    
    navSolutionsKF.usrTime(msIndex,:)        = time;  
    navSolutionsKF.usrPos(msIndex,:)      = estusr_kf;
    navSolutionsKF.usrVel(msIndex,:)      = estVel;
    navSolutionsKF.usrPosENU(msIndex,:)   = usrenu(msIndex,:);
    navSolutionsKF.usrVelENU(msIndex,:)   = usrenuvel(msIndex,:);
    navSolutionsKF.usrPosLLH(msIndex,:)   = usrllh(msIndex,:);
    navSolutionsKF.error_state(msIndex,:)       = error_state;
    navSolutionsKF.clkDrift(msIndex)       = clkDrift;
    navSolutionsKF.clkBias(msIndex)       = clkBias_kf;
    navSolutionsKF.isb(msIndex)           = total_state(9);
    navSolutionsKF.newZ(msIndex,:)        = Z;
    navSolutionsKF.meas_inno(msIndex,:)   = (Z - H * error_state);
    navSolutionsKF.predicted_z(msIndex,:) = H * error_state;
    navSolutionsKF.kalman_gain(:,:,msIndex)	 = kalman_gain;
    navSolutionsKF.state_cov(msIndex,:)      = diag(state_cov);
    navSolutionsKF.satEA(msIndex,:)      = el;
    navSolutionsKF.satAZ(msIndex,:)      = az;
   
    % predict postion and clkBias at next epoch   
	total_state  = Transistion_Matrix * total_state;
    estusr_kf = total_state(1:3)';
    clkBias_kf = total_state(7)';
 
    % adaptively update measurement noise covariance R
    if flag_corrCovEst2 == 1 && counterUptR == thresUptR
        tmpR = diag(var(recordR));
        mesurement_noise(1:svlength,1:svlength)=tmpR(1:svlength,1:svlength) ;% ;  
        mesurement_noise(svlength+1:2*svlength,svlength+1:2*svlength)= ...
                                tmpR(svlength+1:2*svlength,svlength+1:2*svlength) ;% 
        
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
        navSolutionsKF.R(counter_r,:) = diag(mesurement_noise);
    end

    % 
    TOW_USR_CT(msIndex) = time;
    if mod(msIndex, track.progressEvery) == 0
        fprintf('CT_KF: index = %6d TOW = %f E = %f N = %f U = %f  B = %f D = %f\n',msIndex, ...
            TOW_USR_CT(msIndex),usrenu(msIndex,1),usrenu(msIndex,2),usrenu(msIndex,3), clkBias_kf,clkDrift); 
    end
 
end % end for msIndex
navSolutionsKF.refPos = cnslxyz;
navSolutionsKF.sv     = sv;           % channel index (NavIC = 100 + PRN)
navSolutionsKF.prBase0 = prBase0;   % pseudorange time origin, shared with trackingUT_POS
navSolutionsCT = navSolutionsKF;

save(['navSolCT_',file.fileName,'.mat'], 'navSolutionsCT','eph','TOW_USR_CT');
save(['tckRstCT_',file.fileName,'.mat'], 'TckResultCT','CN0_CT'); 

