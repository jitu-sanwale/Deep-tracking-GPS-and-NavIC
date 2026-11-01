function [pseudorange, relative_time, base]= pr_est(Acquired, for_prest,signal, TckResult, msIndex)%
%Purpose
%   pseudorange estimation
%Inputs: 
%	eph             - ephemeris 
%	Acquired        - acquisition results
%	for_prest     	- parameters used for pseudorange estimation
%	signal          - parameters related to signals
%	codedelay_pos 	- code delay used for positioning
%	pdi             - predetection integration interval
%	msIndex         - index of millisecond 
%Outputs:
%	pseudorange         - estimated pseudorange
%	relative_time   	- time difference between different channels
%	base                - fixed receiver time origin (samples), reused by
%	                      trackingUT_POS so its pseudoranges share CT's clock
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU


%% 
sv      = Acquired.sv;
c       = 299792458;
bSec    = 75*1e-3; % base second: the distance between sv and usr is between 67 to 86 sec.
fs      = signal.Fs; 
ms      = 1e-3; 
finefreq = 1/fs;

for svindex = 1 : length(Acquired.sv)
    prn = sv(svindex);
    if length(for_prest.sfb1(prn)) > 0
        sfb1            = for_prest.sfb1(prn); % unit: 20ms, the beginning of subframe 1
        nav1            = for_prest.nav1(prn); % unit: ms, the first navigation data point  
        bca             = TckResult(prn).codedelay(msIndex) -1;  % unit: finefreq, beginning of the C/A code 
%         bca             = TckResult(svindex) -1;  % unit: finefreq, beginning of the C/A code
        dat(svindex)    = 20*fs*ms*sfb1 + fs*ms*nav1 + bca; % unit: finefreq.
    end
end % end of svindex

% Keep the reference (receiver time origin) fixed after the first epoch.
% The original code re-used min(dat) every epoch, so the common term of all
% pseudoranges followed the range rate of the reference satellite instead of
% the receiver clock; over long runs this conflicts with the Doppler-derived
% clock drift in the EKF.
persistent base0
if msIndex == 1 || isempty(base0)
    base0 = min(dat);
end
base = base0;
for svindex = 1 : length(Acquired.sv)
    if length(for_prest.sfb1(Acquired.sv(svindex))) > 0
        diff_of_dat(svindex)   = dat(svindex) - base;
        pseudorange(svindex)   = c * (bSec + diff_of_dat(svindex) * finefreq);
        relative_time(svindex) = diff_of_dat(svindex) * finefreq;
    end
end
