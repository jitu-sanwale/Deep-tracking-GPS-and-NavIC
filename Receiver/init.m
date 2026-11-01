%INIT  NavIC/GPS software receiver: classical (CT), vector (VT) and ultratight (UT).
%   Processes IF data sets produced by the simulator (..\Simulator\init.m) or
%   any data folder with the same layout:
%     <dataDir>\IFdata\gps_0IF_Fs5MHz_8bit_IQ_Fighter_<name>.bin     (required)
%     <dataDir>\IFdata\navic_0IF_Fs5MHz_8bit_IQ_Fighter_<name>.bin   (NavIC configurations)
%     <dataDir>\FighterProfiles\Fighter_<name>.csv                   (truth, for the comparison)
%     <dataDir>\FighterProfiles\Fighter_<name>_IMU_200Hz.mat         (IMU, for UT)
%   The GPS file is the timing reference of the receiver and is needed in
%   every configuration, including NavIC alone.
%
%   Prompts for the data folder, the manoeuvres found in it (one, several,
%   or all), the GNSS configuration (GPS alone, NavIC alone, GPS + NavIC)
%   and the receiver modes, then runs the receiver with the truth comparison
%   and writes, according to the selection:
%     every case       <outDir>\results_<dataset>\  plots (positions, errors,
%                      C/N0, clock, UT diagnostics), error CSVs,
%                      truth_comparison.mat, and the case results table
%                      case_table.md / .tex (also printed)
%     several cases    <outDir>\results_FighterSummary*\  cross-manoeuvre
%                      summary table and bar chart (csv, md, png)
%     all cases        <outDir>\report_all_cases\  when the eight predefined
%                      manoeuvres are available in all three configurations:
%                      result tables for NavIC alone, GPS alone and GPS + NavIC,
%                      profile and satellite-geometry tables, DOPs, summary
%                      statistics, and the comparison and diagnostic figures
%   outDir defaults to Receiver\output.
%
%   Non-interactive use: define a struct RX before running init, e.g.
%     RX.select = {'turn'};  RX.config = {'NavIC'};  init
%   Fields: dataDir, outDir, select (manoeuvre names or 'all'), config
%   (subset of {'GPS','NavIC','GPS+NavIC'}), mode (1 CT, 2 VT, 3 UT, 4 all;
%   default 4), rerun (default false).
%
%   Run time: 17-40 min per 60 s data set and configuration (mode 4).

rxDir = fileparts(mfilename('fullpath'));
addpath(rxDir, fullfile(rxDir, 'acqtckpos'), fullfile(rxDir, 'geo'), fullfile(rxDir, 'common'), fullfile(rxDir, 'report'));
if ~exist('RX', 'var') || ~isstruct(RX), RX = struct(); end
RXRUN = rxConfigure(RX, rxDir);
rxProcess(RXRUN, rxDir);
clear RXRUN

% ======================================================================
function C = rxConfigure(S, rxDir)
CONF = {'GPS', '', 'GPS alone'; 'NavIC', '_NavIC', 'NavIC alone'; 'GPS+NavIC', '_GPS_NavIC', 'GPS + NavIC'};
inter = ~isfield(S, 'select');
C.dataDir = getf(S, 'dataDir', fullfile(fileparts(rxDir), 'Simulator'));
C.outDir  = getf(S, 'outDir', fullfile(rxDir, 'output'));
fprintf('\n=== NavIC / GPS receiver: classical, vector and ultratight tracking ===\n');
if inter
    d = strtrim(strrep(input(sprintf('Data folder (with IFdata and FighterProfiles) [%s]: ', C.dataDir), 's'), '"', ''));
    if ~isempty(d), C.dataDir = d; end
end
[names, hasN, hasT, hasI] = scanData(C.dataDir);
if isempty(names)
    error('init: no data sets (IFdata\\gps_0IF_Fs5MHz_8bit_IQ_Fighter_<name>.bin) in %s - run the simulator first.', C.dataDir);
end
predefined = fighterProfileList();
if inter
    fprintf('\nData sets in %s:\n', C.dataDir);
    fprintf('  %3s  %-22s %-5s %-6s %-6s %s\n', '', 'manoeuvre', 'GPS', 'NavIC', 'truth', 'IMU');
    yn = {'-', 'yes'};
    for k = 1:numel(names)
        tag = '';  if any(strcmp(names{k}, predefined)), tag = '(predefined)'; end
        fprintf('  %3d) %-22s %-5s %-6s %-6s %-4s %s\n', k, names{k}, 'yes', yn{hasN(k)+1}, yn{hasT(k)+1}, yn{hasI(k)+1}, tag);
    end
    fprintf('  %3d) all of them\n', numel(names) + 1);
    fprintf(['  One number processes a single case; several numbers (e.g. "2 4") a combination;\n' ...
             '  %d all cases (with all three configurations: the full report of tables and figures).\n'], numel(names) + 1);
    sel = askNumbers(sprintf('Select [%d]: ', numel(names) + 1), numel(names) + 1, 1:numel(names) + 1);
    if any(sel == numel(names) + 1), C.select = names; else, C.select = names(sel); end
    c = askNumbers('Configuration: 1) GPS alone  2) NavIC alone  3) GPS + NavIC  4) all three [4]: ', 4, 1:4);
    if c == 4, c = 1:3; end
    C.config = CONF(c, 1)';
    C.mode = askNumbers('Receiver modes: 1) classical  2) vector  3) ultratight  4) all three [4]: ', 4, 1:4);
    d = strtrim(strrep(input(sprintf('Output folder [%s]: ', C.outDir), 's'), '"', ''));
    if ~isempty(d), C.outDir = d; end
    C.rerun = askNumbers('Re-run data sets that already have results? 0) no  1) yes [0]: ', 0, 0:1) == 1;
else
    C.select = getf(S, 'select', 'all');
    if ischar(C.select) && strcmpi(C.select, 'all'), C.select = names; end
    C.select = cellstr(C.select);
    C.config = cellstr(getf(S, 'config', {'GPS', 'NavIC', 'GPS+NavIC'}));
    C.mode = getf(S, 'mode', 4);  C.rerun = getf(S, 'rerun', false);
end
[ok, ic] = ismember(C.config, CONF(:,1));
if ~all(ok), error('init: unknown configuration (use GPS, NavIC, GPS+NavIC).'); end
C.suffix = CONF(ic, 2)';  C.label = CONF(ic, 3)';
% checks: NavIC file for the NavIC configurations, IMU for UT, truth for the comparison
for k = 1:numel(C.select)
    i = find(strcmp(C.select{k}, names));
    if isempty(i), error('init: no data set "%s" in %s.', C.select{k}, C.dataDir); end
    if any(contains(C.config, 'NavIC')) && ~hasN(i), error('init: %s has no NavIC IF file.', C.select{k}); end
    if any(C.mode == [3 4]) && ~hasI(i), error('init: %s has no IMU file (needed for ultratight).', C.select{k}); end
    if ~hasT(i), fprintf(2, 'warning: %s has no truth trajectory - no truth comparison.\n', C.select{k}); end
end
C.requested = C.select;   % what the user asked for; reports follow this selection
% warm start: the predefined and the custom manoeuvres of the simulator
% share the baseline's first 12 s, so their runs reuse its decoded ephemeris;
% process the baseline first when it is available
C.warm = false(1, numel(C.select));
if any(strcmp(names, 'baseline_cruise'))
    for k = 1:numel(C.select)
        C.warm(k) = ~strcmp(C.select{k}, 'baseline_cruise') && sharesBaselineStart(C.dataDir, C.select{k});
    end
    if any(C.warm) && ~any(strcmp(C.select, 'baseline_cruise'))
        C.select = [{'baseline_cruise'} C.select];  C.warm = [false C.warm];
        fprintf('baseline_cruise added: its decoded ephemeris is reused by the other runs (warm start).\n');
    end
end
[~, o] = sort(~strcmp(C.select, 'baseline_cruise'));  C.select = C.select(o);  C.warm = C.warm(o);
fprintf('\nPlan: %d data set(s) x %d configuration(s), mode %d -> %d run(s), roughly %d-%d min\n', ...
    numel(C.select), numel(C.config), C.mode, numel(C.select)*numel(C.config), ...
    17*numel(C.select)*numel(C.config), 40*numel(C.select)*numel(C.config));
fprintf('  data sets: %s\n  configurations: %s\n  output: %s\n', strjoin(C.select, ', '), strjoin(C.label, ', '), C.outDir);
if inter && askNumbers('Start? 1) yes  0) no [1]: ', 1, 0:1) == 0, error('init: cancelled.'); end
end

% ----------------------------------------------------------------------
function rxProcess(C, rxDir)
if ~exist(C.outDir, 'dir'), mkdir(C.outDir); end
here = pwd;  back = onCleanup(@() cd(here));
tAll = tic;
for c = 1:numel(C.config)
    for k = 1:numel(C.select)
        p = C.select{k};
        tag = ['gps_0IF_Fs5MHz_8bit_IQ_Fighter_' p C.suffix{c}];
        if exist(fullfile(C.outDir, ['results_' tag], 'truth_comparison.mat'), 'file') && ~C.rerun
            fprintf('[receiver] %s, %s: results exist - kept\n', C.label{c}, p);  continue;
        end
        if C.warm(k)
            ref = ['gps_0IF_Fs5MHz_8bit_IQ_Fighter_baseline_cruise' C.suffix{c}];
            for pre = {'eph', 'sbf', 'TckResultCT_forEph'}
                src = fullfile(C.outDir, [pre{1} '_' ref '.mat']);  dst = fullfile(C.outDir, [pre{1} '_' tag '.mat']);
                if exist(src, 'file') && ~exist(dst, 'file'), copyfile(src, dst); end
            end
        end
        for pre = {'navSolCT', 'navSolVT', 'navSolUT', 'tckRstCT', 'tckRstVT', 'tckRstUT'}
            f = fullfile(C.outDir, [pre{1} '_' tag '.mat']);  if exist(f, 'file'), delete(f); end
        end
        fprintf('\n[receiver] %s, %s\n', C.label{c}, p);
        setenv('VT_RXMODE', num2str(C.mode));
        setenv('VT_GNSS', C.config{c});
        setenv('VT_OUTDIR', C.outDir);
        setenv('VT_DATADIR', C.dataDir);
        setenv('VT_DATAFILE',       fullfile(C.dataDir, 'IFdata', ['gps_0IF_Fs5MHz_8bit_IQ_Fighter_' p '.bin']));
        setenv('VT_DATAFILE_NAVIC', fullfile(C.dataDir, 'IFdata', ['navic_0IF_Fs5MHz_8bit_IQ_Fighter_' p '.bin']));
        setenv('VT_IMUFILE',        fullfile(C.dataDir, 'FighterProfiles', ['Fighter_' p '_IMU_200Hz.mat']));
        setenv('VT_TRUTHFILE',      fullfile(C.dataDir, 'FighterProfiles', ['Fighter_' p '.csv']));
        diary(fullfile(C.outDir, ['run_log_' tag '.txt']));  diary on;
        try
            runSDR(rxDir);
        catch err
            fprintf(2, 'receiver run %s failed: %s\n', tag, err.message);
        end
        diary off;  fclose('all');
    end
end

% results: a table for every case, a cross-manoeuvre summary for several
% cases, and the full report when all cases are available
predefined = fighterProfileList();
fullReport = all(ismember(predefined, C.requested));
for c = 1:numel(C.config)
    for k = 1:numel(C.requested)
        res = fullfile(C.outDir, ['results_gps_0IF_Fs5MHz_8bit_IQ_Fighter_' C.requested{k} C.suffix{c}]);
        caseReport(res, sprintf('%s, %s', C.requested{k}, C.label{c}));
    end
    if numel(C.requested) > 1 && ~fullReport
        S = summarizeFighterRuns(C.requested, C.suffix{c}, C.outDir, C.dataDir); %#ok<NASGU>
    end
end
if fullReport
    sfx = {'', '_NavIC', '_GPS_NavIC'};
    have = false(numel(predefined), numel(sfx));
    for s = 1:numel(sfx)
        for k = 1:numel(predefined)
            have(k, s) = exist(fullfile(C.outDir, ['results_gps_0IF_Fs5MHz_8bit_IQ_Fighter_' predefined{k} sfx{s}], ...
                                        'truth_comparison.mat'), 'file') == 2;
        end
    end
    if all(have(:))
        for s = 1:numel(sfx), summarizeFighterRuns(predefined, sfx{s}, C.outDir, C.dataDir); end
        rep = fullfile(C.outDir, 'report_all_cases');
        reportGeometry(C.outDir, C.dataDir, rep);            % NavIC satellite table and DOPs
        reportProfileTable(predefined, C.dataDir, rep);      % manoeuvre dynamics table
        reportTables(predefined, rep, C.outDir);             % result tables and statistics
        reportFigures(predefined, C.dataDir, rep, C.outDir); % comparison and diagnostic figures
        fprintf('\nFull report (all cases) written to %s\n', rep);
    else
        for c = 1:numel(C.config), summarizeFighterRuns(C.requested, C.suffix{c}, C.outDir, C.dataDir); end
        conf = {'GPS alone', 'NavIC alone', 'GPS + NavIC'};
        fprintf(['\nFull report not written: it needs all eight predefined manoeuvres in all three ' ...
                 'configurations; missing: %s\n'], strjoin(conf(~all(have, 1)), ', '));
    end
end
fprintf('\nReceiver finished in %.1f min. Results: %s\n', toc(tAll)/60, C.outDir);
end

function runSDR(rxDir)
% SDR_main is a script that clears its workspace; running it inside this
% function keeps the clear local to the call
run(fullfile(rxDir, 'SDR_main.m'));
end

% ----------------------------------------------------------------------
function [names, hasN, hasT, hasI] = scanData(dataDir)
d = dir(fullfile(dataDir, 'IFdata', 'gps_0IF_Fs5MHz_8bit_IQ_Fighter_*.bin'));
names = regexprep({d.name}, '^gps_0IF_Fs5MHz_8bit_IQ_Fighter_(.*)\.bin$', '$1');
predefined = fighterProfileList();
[~, o] = sort(cellfun(@(n) find([strcmp(n, predefined) true], 1), names));  names = names(o);   % predefined order first
f = @(sub, pat) cellfun(@(n) exist(fullfile(dataDir, sub, sprintf(pat, n)), 'file') > 0, names);
hasN = f('IFdata', 'navic_0IF_Fs5MHz_8bit_IQ_Fighter_%s.bin');
hasT = f('FighterProfiles', 'Fighter_%s.csv');
hasI = f('FighterProfiles', 'Fighter_%s_IMU_200Hz.mat');
end

function tf = sharesBaselineStart(dataDir, name)
% predefined and custom-builder manoeuvres (6-DOF model) start like the baseline;
% imported CSV trajectories (model 'trajectory') may not
tf = false;
f = fullfile(dataDir, 'FighterProfiles', ['Fighter_' name '_IMU_200Hz.mat']);
if ~exist(f, 'file'), return; end
S = load(f, 'imu');
tf = isfield(S.imu, 'profile') && isfield(S.imu.profile, 'stats') && isfield(S.imu.profile.stats, 'model') ...
     && strcmp(S.imu.profile.stats.model, '6dof');
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
