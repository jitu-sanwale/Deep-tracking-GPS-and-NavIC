%INIT  NavIC/GPS IF signal simulator for high-performance aircraft (standalone).
%   Generates, for the manoeuvres you select,
%     FighterProfiles\Fighter_<name>.csv            1 kHz truth trajectory
%     FighterProfiles\Fighter_<name>_IMU_200Hz.*    navigation-grade IMU data + truth
%     IFdata\gps_0IF_Fs5MHz_8bit_IQ_Fighter_<name>.bin    GPS L1 C/A
%     IFdata\navic_0IF_Fs5MHz_8bit_IQ_Fighter_<name>.bin  NavIC L5 SPS (+ _time.mat)
%   in the data folder (default: this folder). IF files: 5 MHz complex
%   baseband, 8-bit I/Q, 10 MB per second of data (600 MB per 60 s file).
%
%   Choices (prompted):
%     1-8  predefined manoeuvres (F-16 6-DOF model + NDI autopilot)
%     9    all eight predefined manoeuvres
%     10   custom manoeuvre, assembled from building blocks (turn, bank,
%          roll, pull/push, loop, ...) and flown by the same 6-DOF model
%     11   custom trajectory from your own CSV file
%   Several numbers may be entered at once, e.g. "2 4 10".
%
%   Non-interactive use: define a struct SIM before running init, e.g.
%     SIM.select = {'turn','break_gonset'};  SIM.gnss = 'both';  init
%     SIM.select = {'custom_myturn'};  SIM.custom.custom_myturn = {{'turn',6,8,3,1},{'roll',360,270}};  init
%     SIM.select = {'custom_track'};   SIM.csv.custom_track = 'C:\data\track.csv';  init
%   Fields: select, gnss ('GPS'|'NavIC'|'both'), seconds (default 60),
%   overwrite (default false), dataDir, custom.<name> (block list),
%   csv.<name> (CSV path), description.<name>.
%
%   Defaults: 60 s, both bands, IMU seed 1, +/-20 g accelerometer range.
%   The receiver (..\Receiver\init.m) reads the data folder written here.

simDir = fileparts(mfilename('fullpath'));
addpath(simDir, fullfile(simDir, 'F16Aircraft'));
if ~exist('SIM', 'var') || ~isstruct(SIM), SIM = struct(); end
SIMRUN = simConfigure(SIM, simDir);
simGenerate(SIMRUN, simDir);
SIMDATADIR = SIMRUN.dataDir;   % for ..\init.m (generate, then process)
clear SIMRUN

% ======================================================================
function C = simConfigure(S, simDir)
predefined = fighterProfileList();
inter = ~isfield(S, 'select');
C.dataDir = getf(S, 'dataDir', simDir);
C.custom = getf(S, 'custom', struct());  C.csv = getf(S, 'csv', struct());
C.description = getf(S, 'description', struct());

fprintf('\n=== NavIC / GPS IF simulator for high-performance aircraft ===\n');
if inter
    fprintf('\nPredefined manoeuvres:\n');
    for k = 1:numel(predefined)
        P = fighterPhases(predefined{k});
        fprintf('  %2d) %-17s %s\n', k, predefined{k}, P.description);
    end
    fprintf('   9) all eight predefined manoeuvres\n');
    fprintf('  10) custom manoeuvre (building blocks, 6-DOF model)\n');
    fprintf('  11) custom trajectory from a CSV file\n');
    sel = askNumbers('Select (e.g. "9", or "2 4 10") [9]: ', 9, 1:11);
    C.select = {};
    for s = sel(:)'
        if s <= 8,      C.select{end+1} = predefined{s};
        elseif s == 9,  C.select = [C.select predefined];
        elseif s == 10
            [nm, blocks, desc] = buildCustomManeuver();
            C.select{end+1} = nm;  C.custom.(nm) = blocks;  C.description.(nm) = desc;
        else
            [nm, f] = askCsv();
            C.select{end+1} = nm;  C.csv.(nm) = f;
        end
    end
    C.select = unique(C.select, 'stable');
    b = askNumbers('GNSS band: 1) GPS L1  2) NavIC L5  3) both (needed for every receiver configuration) [3]: ', 3, 1:3);
    C.gnss = {'GPS', 'NavIC', 'both'};  C.gnss = C.gnss{b};
    C.seconds = askNumbers('Record length in seconds (13-120) [60]: ', 60, []);
    d = input(sprintf('Data folder [%s]: ', C.dataDir), 's');
    if ~isempty(strtrim(d)), C.dataDir = strtrim(d); end
    C.overwrite = askNumbers('Overwrite files that already exist? 0) no, skip them  1) yes [0]: ', 0, 0:1) == 1;
else
    C.select = cellstr(S.select);
    C.gnss = getf(S, 'gnss', 'both');
    C.seconds = getf(S, 'seconds', 60);
    C.overwrite = getf(S, 'overwrite', false);
end
C.seconds = min(max(C.seconds, 13), 120);
for k = 1:numel(C.select)
    nm = C.select{k};
    if ~any(strcmp(nm, predefined)) && ~isfield(C.custom, nm) && ~isfield(C.csv, nm)
        error('init: "%s" is neither a predefined manoeuvre nor defined in SIM.custom / SIM.csv.', nm);
    end
end
if strcmpi(C.gnss, 'both'), C.bands = {'GPS', 'NavIC'}; else, C.bands = {C.gnss}; end
nF = numel(C.select)*numel(C.bands);
fprintf('\nPlan: %d manoeuvre(s) x %s, %g s -> %d IF file(s), about %.1f GB, roughly %d-%d min\n', ...
    numel(C.select), strjoin(C.bands, ' + '), C.seconds, nF, nF*C.seconds*0.01, ...
    round(nF*C.seconds/60*20), round(nF*C.seconds/60*30));
fprintf('  manoeuvres: %s\n  data folder: %s\n', strjoin(C.select, ', '), C.dataDir);
if inter
    if askNumbers('Start? 1) yes  0) no [1]: ', 1, 0:1) == 0, error('init: cancelled.'); end
end
end

% ----------------------------------------------------------------------
function simGenerate(C, simDir)
profDir = fullfile(C.dataDir, 'FighterProfiles');  ifDir = fullfile(C.dataDir, 'IFdata');
if ~exist(profDir, 'dir'), mkdir(profDir); end
if ~exist(ifDir, 'dir'), mkdir(ifDir); end
predefined = fighterProfileList();
tAll = tic;
for k = 1:numel(C.select)
    nm = C.select{k};
    base = fullfile(profDir, ['Fighter_' nm]);
    have = exist([base '.csv'], 'file') && exist([base '_IMU_200Hz.mat'], 'file');
    if ~have || C.overwrite
        fprintf('\n[profile] %s\n', nm);
        if isfield(C.csv, nm)
            importCustomTrajectory(C.csv.(nm), nm, struct('outDir', profDir));
        else
            opt = struct('outDir', profDir, 'accelRangeG', 20, 'seed', 1, 'T', max(C.seconds, 60));
            if ~any(strcmp(nm, predefined))
                opt.customPhases = C.custom.(nm);
                if isfield(C.description, nm), opt.customDescription = C.description.(nm); end
                blocks = C.custom.(nm);  save([base '_spec.mat'], 'blocks');   % record of the definition
            end
            generateFighterProfile(nm, opt);
        end
    else
        fprintf('\n[profile] %s exists - kept\n', nm);
    end
    for b = 1:numel(C.bands)
        f = fullfile(ifDir, sprintf('%s_0IF_Fs5MHz_8bit_IQ_Fighter_%s.bin', lower(C.bands{b}), nm));
        d = dir(f);
        if ~isempty(d) && ~C.overwrite && d.bytes >= C.seconds*1e7
            fprintf('[IF] %s exists - kept\n', f);  continue;
        end
        if ~isempty(d), delete(f); end          % never write through an old file (may be a hard link)
        generateIFdata(nm, C.bands{b}, struct('profileDir', profDir, 'ifDir', ifDir, 'ms', round(C.seconds*1000)));
    end
end
fprintf('\nSimulator finished in %.1f min. Data folder: %s\n', toc(tAll)/60, C.dataDir);
end

% ----------------------------------------------------------------------
function [nm, blocks, desc] = buildCustomManeuver()
B = customManeuverBlocks();
fprintf('\nCustom manoeuvre. The record starts with 12 s of straight and level flight at\n');
fprintf('250 m/s, 5000 m, heading East (as all predefined manoeuvres); the blocks follow in order,\n');
fprintf('then the aircraft levels off and cruises to the end of the record.\n');
showBlocks(B);
fprintf('Enter one block per line, e.g. "turn 7 10 3 1"; omitted parameters take the default.\n');
fprintf('Commands: "list" shows the blocks, "undo" removes the last one, empty line finishes.\n');
blocks = {};
while true
    s = strtrim(input(sprintf('block %d> ', numel(blocks) + 1), 's'));
    if isempty(s)
        if isempty(blocks), fprintf('  at least one block is needed\n'); continue; end
        break;
    end
    if strcmpi(s, 'list'), showBlocks(B); continue; end
    if strcmpi(s, 'undo'), if ~isempty(blocks), blocks(end) = []; end, continue; end
    tok = strsplit(s);  i = find(strcmpi(tok{1}, {B.name}));
    if isempty(i), fprintf('  unknown block "%s" (type "list")\n', tok{1}); continue; end
    v = str2double(tok(2:end));
    if any(isnan(v)) || numel(v) > numel(B(i).params)
        fprintf('  %s takes up to %d numbers: %s\n', B(i).name, numel(B(i).params), strjoin(B(i).params, ', ')); continue;
    end
    a = B(i).defaults;  a(1:numel(v)) = v;
    bad = a < B(i).lo | a > B(i).hi;
    if any(bad)
        j = find(bad, 1);
        fprintf('  %s must lie in [%g, %g]\n', B(i).params{j}, B(i).lo(j), B(i).hi(j)); continue;
    end
    blocks{end+1} = [{B(i).name}, num2cell(a)]; %#ok<AGROW>
    fprintf('  added: %s %s\n', B(i).name, mat2str(a));
end
nm = strtrim(input('Name for this manoeuvre [custom_1]: ', 's'));
if isempty(nm), nm = 'custom_1'; end
nm = cleanName(nm);
desc = strtrim(input('Short description (optional): ', 's'));
if isempty(desc), desc = strjoin(cellfun(@(b) sprintf('%s %s', b{1}, mat2str([b{2:end}])), blocks, 'UniformOutput', false), '; '); end
fprintf('  saved as %s. To repeat it without prompts:\n  SIM.custom.%s = %s;\n', nm, nm, blocksToString(blocks));
end

function showBlocks(B)
fprintf('\n  %-7s %-40s %s\n', 'block', 'parameters (defaults)', 'meaning');
for i = 1:numel(B)
    p = strjoin(cellfun(@(n, d) sprintf('%s=%g', n, d), B(i).params, num2cell(B(i).defaults), 'UniformOutput', false), ' ');
    fprintf('  %-7s %-40s %s\n', B(i).name, p, B(i).help);
end
fprintf('\n');
end

function s = blocksToString(blocks)
c = cellfun(@(b) ['{''' b{1} ''',' strjoin(cellfun(@(x) sprintf('%g', x), b(2:end), 'UniformOutput', false), ',') '}'], ...
            blocks, 'UniformOutput', false);
s = ['{' strjoin(c, ', ') '}'];
end

function [nm, f] = askCsv()
fprintf('\nCustom trajectory CSV: one header line, then columns\n');
fprintf('  Time [sec], Latitude [deg], Longitude [deg], Height [meter], Roll [deg], Pitch [deg], Heading [deg]\n');
fprintf('  any sample rate, at least 13 s; it is resampled to 1 kHz and an IMU record is synthesised.\n');
while true
    f = strtrim(strrep(input('CSV file: ', 's'), '"', ''));
    if exist(f, 'file'), break; end
    fprintf('  file not found\n');
end
[~, b] = fileparts(f);
nm = strtrim(input(sprintf('Name for this trajectory [custom_%s]: ', cleanName(b, false)), 's'));
if isempty(nm), nm = ['custom_' cleanName(b, false)]; end
nm = cleanName(nm);
end

function nm = cleanName(nm, prefix)
if nargin < 2, prefix = true; end
nm = regexprep(nm, '[^A-Za-z0-9_]', '_');
if prefix && ~strncmp(nm, 'custom', 6), nm = ['custom_' nm]; end
end

function v = askNumbers(prompt, def, allowed)
while true
    s = strtrim(input(prompt, 's'));
    if isempty(s), v = def; return; end
    v = str2double(strsplit(s));
    if all(~isnan(v)) && (isempty(allowed) || all(ismember(v, allowed))), return; end
    fprintf('  please enter valid number(s)\n');
end
end

function v = getf(S, f, def)
if isfield(S, f) && ~isempty(S.(f)), v = S.(f); else, v = def; end
end
