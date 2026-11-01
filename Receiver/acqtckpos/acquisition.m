function  Acquired = acquisition(file,signal,acq)
%Purpose:
%   Perform signal acquisition using parallel code phase search algorithm
%   Fine frequency is found using long FFT.
%Inputs:
%	file        - parameters related to the data file to be processed,a structure
%	signal      - parameters related to signals,a structure
%Outputs:
%	Acquired	- acquisition results, i.e. space vehicle (SV),signal-to-noise 
%                   ratio (SNR), Doppler, code delay, and fine carrier
%                   frequency
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU
% Modified: generic real / complex (I/Q) data reader, complex carrier
% wipe-off and signed fine-frequency search for zero-IF complex data.


%% Initialization 
Acquired.sv = []; 
Acquired.SNR = []; 
Acquired.Doppler = []; 
Acquired.codedelay = []; 
Acquired.fineFreq = [];
sampleindex = 1: signal.Sample;
mixSign = 1;
if isfield(signal,'mixSign'), mixSign = signal.mixSign; end
bytesPerSample = file.dataPrecision*file.dataType;

%% read data (GPS L1 file, and the NavIC L5 file when it is used)
fseek(file.fid,file.skip*signal.Sample*bytesPerSample,'bof');
rawsignalG = readIFData(file.fid, signal.Sample*4, file); % 4ms data for acquisition
chanList = acq.prnList;   % empty for a NavIC-only receiver
if isfield(file,'fidN') && file.fidN > 0
    fseek(file.fidN,file.skip*signal.Sample*bytesPerSample,'bof');
    rawsignalN = readIFData(file.fidN, signal.Sample*4, file);
    chanList = [chanList, 100 + acq.prnListN];   % NavIC channel index = 100 + PRN
end

%% acquisition process 
carrier = zeros(acq.freqNum, signal.Sample);
for freqband = 1 : acq.freqNum
    dopplershift = acq.freqMin + acq.freqStep*(freqband-1);
    if file.dataType == 2
        carrier(freqband,:) = exp(1i*mixSign*2*pi*(signal.IF + dopplershift) * sampleindex ./ signal.Fs);
    else
        carrier(freqband,:) = exp(1i*2*pi*(signal.IF + dopplershift) * sampleindex ./ signal.Fs); 
    end
end

for svindex = chanList
    ch = gnssChan(svindex, signal);
    ocode = ch.code;
    if ch.isNavIC, rawsignal = rawsignalN; else, rawsignal = rawsignalG; end
    ocode = [ocode ocode];
    scode = ocode(ceil(sampleindex.*(signal.codeFreqBasis/signal.Fs)));
    temp3 = fft(scode);
    correlation = zeros(acq.freqNum,signal.Sample);
    for idx = 1 : 4  % 4ms data
        blk = rawsignal(1+(idx-1)*signal.Sample:idx*signal.Sample);
        for freqband = 1 : acq.freqNum        
            temp1 = blk.* carrier(freqband,:);
            temp2 = conj(fft(temp1));
            correlation(freqband,:) = correlation(freqband,:) + abs(ifft(temp3.*temp2)).^2;
        end
    end
    [peak, fbin] = max(max(abs(correlation')));
    [peak, codePhase] = max(max(abs(correlation)));
    Doppler = acq.freqMin + acq.freqStep * (fbin-1);
    
    codechipshift = ceil(signal.Fs/signal.codeFreqBasis);
    noiseIdx = [1:codePhase-codechipshift codePhase+codechipshift:signal.Sample];
    noiseIdx = noiseIdx(noiseIdx>=1 & noiseIdx<=signal.Sample);
    SNR = 10 * log10((peak^2)/(sum(correlation(fbin,noiseIdx).^2)/length(noiseIdx)));

    if SNR >= acq.threshold  % acquisition threshold
        Acquired.sv        = [Acquired.sv        svindex];
        Acquired.SNR       = [Acquired.SNR       SNR];
        Acquired.Doppler   = [Acquired.Doppler   Doppler];
        Acquired.codedelay = [Acquired.codedelay codePhase-1];
    end
end


%% fine frequency calculation 
fseek(file.fid,file.skip*signal.Sample*bytesPerSample,'bof');
longrawsignalG = readIFData(file.fid, signal.Sample*(acq.L+1), file);
if isfield(file,'fidN') && file.fidN > 0
    fseek(file.fidN,file.skip*signal.Sample*bytesPerSample,'bof');
    longrawsignalN = readIFData(file.fidN, signal.Sample*(acq.L+1), file);
end

for svindex = 1 : length(Acquired.sv)     
    ch = gnssChan(Acquired.sv(svindex), signal);
    caCode = ch.code;
    if ch.isNavIC, longrawsignal = longrawsignalN; else, longrawsignal = longrawsignalG; end
    codeValueIndex = floor((1/signal.Fs*(1:acq.L*signal.Sample))/(1/signal.codeFreqBasis));    
    longCaCode = caCode((rem(codeValueIndex, signal.codelength) + 1));
    CarrSignal = longrawsignal((signal.Sample-Acquired.codedelay(svindex)):...
                    (signal.Sample-Acquired.codedelay(svindex))+acq.L*signal.Sample - 1).* longCaCode;
     
    fftlength = length(CarrSignal) * 20;
    fftSignal = abs(fft(CarrSignal, fftlength));
    
    if file.dataType == 2
        % complex data: search the whole spectrum, keep the sign of the frequency
        [~, FreqPeakIndex] = max(fftSignal);
        f = (FreqPeakIndex-1) * (signal.Fs/fftlength);
        if f >= signal.Fs/2, f = f - signal.Fs; end
        % tracking replica frequency (= IF + Doppler for mixSign = -1)
        fineDoppler = -mixSign * f;
        % sanity check against the coarse search
        if abs((fineDoppler - signal.IF) - Acquired.Doppler(svindex)) > 2*acq.freqStep
            fineDoppler = signal.IF + Acquired.Doppler(svindex);
        end
    else
        halffftlength = ceil((fftlength)/2);
        [~, FreqPeakIndex] = max(fftSignal(1:halffftlength));
        fineDoppler = FreqPeakIndex * (signal.Fs/fftlength);    
    end
    
    Acquired.fineFreq = [Acquired.fineFreq  fineDoppler];
        
    fprintf('SV[%2d] SNR = %2.2f Code phase = %5d  Doppler frequency = %5f\n', ...
            Acquired.sv(svindex),Acquired.SNR(svindex),Acquired.codedelay(svindex), ...
            Acquired.fineFreq(svindex)-signal.IF);
end

end % end function
