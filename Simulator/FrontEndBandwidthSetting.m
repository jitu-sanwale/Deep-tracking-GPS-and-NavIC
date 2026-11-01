function feBandwidth = FrontEndBandwidthSetting(settings, defaultBandwidth)
% -------------------------------------------------------------------------
%                   SoftSim: GPS IF signal simulator
% -------------------------------------------------------------------------
% Prompts for the front-end (FE) filter bandwidth [Hz], constrained to the
% range the current signal configuration can actually support, and
% suggests a sensible default.
%
%   feBandwidth = FrontEndBandwidthSetting(settings, defaultBandwidth)
%
%   Range:
%     - Minimum: 2 x codeFreqBasis, i.e. the null-to-null width of the
%       GPS C/A code's main spectral lobe. Filtering narrower than this
%       clips the main lobe and badly degrades correlation SNR.
%     - Maximum: the widest bandwidth the current samplingFreq/IF/fileType
%       combination can represent without handing getFliterCoef.m/fir1 an
%       invalid (out of (0,1)) normalized cutoff.
%
%   Optional 2nd arg (settings.defaultFeBandwidth) skips the interactive
%   prompt for scripted/batch runs. Leave empty ([]) to always prompt.

minBW = 2 * settings.codeFreqBasis;

if settings.fileType == 1
    % Real IF samples: passband [IF-BW/2, IF+BW/2] must fit within
    % [0, samplingFreq/2].
    maxBW = 2 * min(settings.IF, settings.samplingFreq/2 - settings.IF);
else
    % Complex/IQ baseband samples: usable bandwidth is the full Nyquist
    % span, +/- samplingFreq/2 around 0 Hz.
    maxBW = settings.samplingFreq;
end

if maxBW <= minBW
    error(['FrontEndBandwidthSetting: no valid front-end bandwidth exists ' ...
        'for the current samplingFreq/IF/codeFreqBasis (need max %.3f MHz ' ...
        '> min %.3f MHz). Increase samplingFreq or check settings.IF.'], ...
        maxBW/1e6, minBW/1e6);
end

% Suggested default: the original tool's own choice (2.6 MHz) if it fits
% the current bounds, otherwise the midpoint of the valid range.
suggested = 2.6e6;
if suggested < minBW || suggested > maxBW
    suggested = (minBW + maxBW) / 2;
end

if nargin > 1 && ~isempty(defaultBandwidth)
    feBandwidth = defaultBandwidth;
    if feBandwidth < minBW || feBandwidth > maxBW
        error(['FrontEndBandwidthSetting: settings.defaultFeBandwidth = %.3f MHz ' ...
            'is outside the valid range [%.3f, %.3f] MHz.'], ...
            feBandwidth/1e6, minBW/1e6, maxBW/1e6);
    end
else
    fprintf('\nFront-end filter bandwidth must be between %.3f MHz (GPS C/A main lobe) ', minBW/1e6);
    fprintf('and %.3f MHz (Nyquist limit for the current sampling rate/IF).\n', maxBW/1e6);
    feBandwidth = input(sprintf('Enter front-end filter bandwidth in Hz [default %.3f MHz, press Enter to accept]: ', suggested/1e6));
    if isempty(feBandwidth)
        feBandwidth = suggested;
    end
    while ~(isscalar(feBandwidth) && feBandwidth >= minBW && feBandwidth <= maxBW)
        feBandwidth = input(sprintf('Please enter a value between %.0f and %.0f Hz (or press Enter for %.0f Hz): ', ...
            minBW, maxBW, suggested));
        if isempty(feBandwidth)
            feBandwidth = suggested;
        end
    end
end
