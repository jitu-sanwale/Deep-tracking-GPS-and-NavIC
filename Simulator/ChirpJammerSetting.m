function [chirpEnable, jnrDB, sweepPeriodMs, deltaF0, deltaF1, injStartSec, injEndSec] = ChirpJammerSetting(...
    settings, defaultEnable, defaultSweepPeriodMs, defaultInjStartSec, defaultInjEndSec)
% -------------------------------------------------------------------------
%                   SoftSim: GPS IF signal simulator
% -------------------------------------------------------------------------
% Prompts for a linear-chirp ("sweep") jamming interference, modeled per:
%   Sun, K.; Elhajj, M.; Ochieng, W.Y. "A GNSS Anti-Interference Method
%   Based on Fractional Fourier Transform." IEEE Trans. Aerosp. Electron.
%   Syst. 2024, 60(5), 5636-5651. doi:10.1109/TAES.2024.3366147
%
%   j_RF(t)  = A_inst(t) * cos[2*pi*integral(finst(tau),0,t) + phi0]
%   finst(t) = f0 + k*t,  0 <= t <= tj      (resets every sweep period tj)
%   f0 = fL1 + deltaF0,  f1 = fL1 + deltaF1,  k = (f1-f0)/tj
%
% JNR (jammer-to-noise ratio, dB) and the sweep frequency deviation
% (deltaF0/deltaF1) are NOT prompted here - they are fixed configuration
% read directly from settings.chirpJnrDB/settings.chirpDeltaF0/
% settings.chirpDeltaF1 (set in initSettings.m).
%
% The jammer can also be confined to a start/end time window (elapsed
% seconds from the start of the simulation run, i.e. relative to
% trajectory/output-file time 0, not absolute GPS time) rather than
% running for the whole pass.
%
% Optional args (settings.defaultChirp*) skip the interactive prompts for
% scripted/batch runs. Leave empty ([]) to always prompt.

if nargin > 1 && ~isempty(defaultEnable)
    chirpEnable = defaultEnable;
else
    fprintf('\nChirp (sweep) jamming interference (Sun, Elhajj & Ochieng, IEEE TAES 2024):\n');
    chirpEnable = input('Enable chirp jammer? Enter 1 to enable, 0 to disable: ');
    while ~(isscalar(chirpEnable) && ismember(chirpEnable, [0 1]))
        chirpEnable = input('Please enter 1 (enable) or 0 (disable): ');
    end
end

jnrDB = []; sweepPeriodMs = []; deltaF0 = []; deltaF1 = [];
injStartSec = []; injEndSec = [];
if ~chirpEnable
    return
end

% JNR and sweep frequency deviation come straight from settings (fixed,
% not prompted). Still validate the deviation against the Nyquist limit
% for the current sampling rate, since these are user-edited values.
jnrDB   = settings.chirpJnrDB;
deltaF0 = settings.chirpDeltaF0;
deltaF1 = settings.chirpDeltaF1;
maxDevHz = settings.samplingFreq/2;
if ~(abs(deltaF0) <= maxDevHz && abs(deltaF1) <= maxDevHz)
    error(['ChirpJammerSetting: settings.chirpDeltaF0/chirpDeltaF1 must be ' ...
        'within +/- %.0f Hz (Nyquist) of L1 center for the current ' ...
        'settings.samplingFreq.'], maxDevHz);
end

if nargin > 2 && ~isempty(defaultSweepPeriodMs)
    sweepPeriodMs = defaultSweepPeriodMs;
else
    sweepPeriodMs = input('Enter sweep period in ms [paper used 0.02-0.5 ms; press Enter for default 0.5 ms]: ');
    if isempty(sweepPeriodMs)
        sweepPeriodMs = 0.5;
    end
    while ~(isscalar(sweepPeriodMs) && sweepPeriodMs > 0)
        sweepPeriodMs = input('Please enter a positive sweep period in ms: ');
    end
end

% Injection time window (elapsed seconds from the start of the
% simulation run). Defaults to the whole pass.
totalDurationSec = settings.msToProcess/1000;

if nargin > 3 && ~isempty(defaultInjStartSec)
    injStartSec = defaultInjStartSec;
else
    fprintf('Simulation duration is %.3f s.\n', totalDurationSec);
    injStartSec = input(sprintf('Enter interference start time in s [press Enter for default 0]: '));
    if isempty(injStartSec)
        injStartSec = 0;
    end
    while ~(isscalar(injStartSec) && injStartSec >= 0 && injStartSec < totalDurationSec)
        injStartSec = input(sprintf('Please enter a start time in [0, %.3f) s: ', totalDurationSec));
    end
end

if nargin > 4 && ~isempty(defaultInjEndSec)
    injEndSec = defaultInjEndSec;
else
    injEndSec = input(sprintf('Enter interference end time in s [press Enter for default %.3f, i.e. end of run]: ', totalDurationSec));
    if isempty(injEndSec)
        injEndSec = totalDurationSec;
    end
    while ~(isscalar(injEndSec) && injEndSec > injStartSec && injEndSec <= totalDurationSec)
        injEndSec = input(sprintf('Please enter an end time in (%.3f, %.3f] s: ', injStartSec, totalDurationSec));
    end
end
