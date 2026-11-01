function ifFile = generateIFdata(name, gnss, opt)
%GENERATEIFDATA IF baseband file for one trajectory and one GNSS band.
%   ifFile = generateIFdata(name, gnss, opt)
%
%   name - profile name; the trajectory is <opt.profileDir>\Fighter_<name>.csv
%          (1 kHz, made by generateFighterProfile.m or importCustomTrajectory.m)
%   gnss - 'GPS' (L1 C/A) or 'NavIC' (L5 SPS)
%   opt  - optional struct:
%            profileDir  trajectory folder   (default <simulator>\FighterProfiles)
%            ifDir       output folder       (default <simulator>\IFdata)
%            ms          record length [ms]  (default: whole trajectory, max 60000)
%            seed        noise seed (default 0, MATLAB's default generator state,
%                        which makes regenerated IF files reproducible byte for byte)
%
%   Writes <ifDir>\gps_0IF_Fs5MHz_8bit_IQ_Fighter_<name>.bin or
%   <ifDir>\navic_0IF_Fs5MHz_8bit_IQ_Fighter_<name>.bin (5 MHz complex,
%   8-bit I/Q), plus <...>_time.mat (NavIC time-model truth) for NavIC.
%   All visible satellites; no multipath, jammer or front-end filter by
%   default. The attitude columns drive the body-fixed Rx antenna gain.
%   The GPS and NavIC files of a profile share start time, sample grid,
%   noise, antenna and atmosphere models, so the receiver can combine them.

simDir = fileparts(mfilename('fullpath'));
if nargin < 3, opt = struct(); end
if ~isfield(opt, 'profileDir'), opt.profileDir = fullfile(simDir, 'FighterProfiles'); end
if ~isfield(opt, 'ifDir'),      opt.ifDir      = fullfile(simDir, 'IFdata'); end
if ~isfield(opt, 'ms'),         opt.ms         = []; end
if ~isfield(opt, 'seed'),       opt.seed       = 0; end
if ~any(strcmpi(gnss, {'GPS', 'NavIC'})), error('generateIFdata: gnss must be ''GPS'' or ''NavIC''.'); end

trajFile = fullfile(opt.profileDir, ['Fighter_' name '.csv']);
if ~exist(trajFile, 'file'), error('generateIFdata: trajectory %s not found.', trajFile); end
if ~exist(opt.ifDir, 'dir'), mkdir(opt.ifDir); end

here = pwd;  back = onCleanup(@() cd(here));
cd(simDir);                                   % initSettings uses files relative to the simulator

settings = initSettings();
settings.defaultTrajFile        = trajFile;
settings.defaultSingleSatEnable = 0;
settings.defaultMultipathEnable = 0;
settings.defaultChirpEnable     = 0;
settings.defaultFilterEnable    = 0;
settings.gnss                   = gnss;
settings.IfFile = fullfile(opt.ifDir, sprintf('%s_0IF_Fs5MHz_8bit_IQ_Fighter_%s.bin', lower(gnss), name));
tr = readmatrix(trajFile, 'NumHeaderLines', 1);
msTraj = floor((tr(end,1) - tr(1,1))*1000);
if isempty(opt.ms), settings.msToProcess = min(60000, msTraj);
else,               settings.msToProcess = min(opt.ms, msTraj);
end

fprintf('=== %s IF for %s, %.1f s -> %s\n', gnss, name, settings.msToProcess/1000, settings.IfFile);
tGen = tic;
rngState = rng;  rng(opt.seed, 'twister');      % receiver-noise realisation
signalGenerating(settings);
rng(rngState);
fprintf('=== done in %.1f min\n', toc(tGen)/60);
ifFile = settings.IfFile;
end
