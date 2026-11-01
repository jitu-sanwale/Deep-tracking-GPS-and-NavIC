function reportFigures(profiles, simDir, outDir, runDir)
%REPORTFIGURES Figures of the full test matrix (all predefined manoeuvres,
%   all three configurations):
%     antenna_attitude.png   antenna pattern and LOS versus attitude
%     profiles_dynamics.png  load factor and roll rate of six maneuvers
%     navic_compare.png      maneuver-window rms errors, all configurations
%     break_ut_diag.png      ultratight filter on NavIC alone, break turn
%     turn_error.png         NavIC alone, 3-D errors in the turn
%   (the last two are redrawn from the saved NavIC-only navigation solutions)
%   simDir is the data folder (IFdata\ and FighterProfiles\); runDir the receiver output.
rxDir  = fileparts(fileparts(mfilename('fullpath')));
if nargin < 4 || isempty(runDir), runDir = fullfile(rxDir, 'output'); end
if ~exist(outDir, 'dir'), mkdir(outDir); end
figAntennaAttitude(fullfile(rxDir, 'common'), outDir);
figProfilesDynamics(simDir, outDir);
figNavicCompare(profiles, runDir, outDir);
figTruth('turn',         'T3_error_magnitude.png',        'turn_error.png',    simDir, runDir, outDir);
figTruth('break_gonset', 'T4_ultratight_diagnostics.png', 'break_ut_diag.png', simDir, runDir, outDir);
end

% ---------------------------------------------------------------- rms comparison
function figNavicCompare(profiles, runDir, outDir)
lab = {'Baseline','Turn','Loop','Break','Jinking','Rolls','Split-S/Imm.','Departure'};
rd = @(s) readtable(fullfile(runDir, ['results_FighterSummary' s], 'fighter_summary.csv'), 'TextType', 'string');
G = rd('');  N = rd('_NavIC');  GN = rd('_GPS_NavIC');
series = {N,'CT'; N,'VT'; N,'UT'; GN,'UT'; G,'UT'};
names  = {'CT, NavIC','VT, NavIC','UT, NavIC','UT, NavIC+GPS','UT, GPS'};
col = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; 0.043 0.420 0.282; 0.725 0.725 0.725];
cols = {'man_pos_rms_m','man_vel_rms_mps'};  yl = {'3-D position rms (m)','3-D velocity rms (m/s)'};
ylims = {[1e-1 1e5], [5e-3 1e4]};
f = figure('Color','w','Units','inches','Position',[0.5 0.5 11 6.4],'Visible','off');
tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for q = 1:2
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    set(ax, 'YScale','log', 'FontSize',10, 'LineWidth',1.4, 'TickDir','out', 'XColor','k', 'YColor','k', ...
        'YMinorGrid','on', 'GridAlpha',0.25, 'MinorGridAlpha',0.12);
    M = nan(numel(profiles), size(series,1));
    for i = 1:numel(profiles)
        for s = 1:size(series,1)
            T = series{s,1};  r = T.profile == profiles{i} & T.mode == series{s,2};
            if any(r), M(i,s) = T.(cols{q})(find(r,1)); end
        end
    end
    M(~isfinite(M) | M > 1e6) = NaN;
    b = bar(ax, M, 'grouped', 'EdgeColor', 'k', 'LineWidth', 0.8);
    for s = 1:numel(b), b(s).FaceColor = col(s,:); end
    set(ax, 'XTick', 1:numel(profiles), 'XTickLabel', lab(1:numel(profiles)));
    ylabel(ax, yl{q}); ylim(ax, ylims{q});
    if q == 1
        legend(ax, names, 'Location', 'northwest', 'Box', 'off', 'NumColumns', 5);
    end
end
exportgraphics(f, fullfile(outDir, 'navic_compare.png'), 'Resolution', 200);  close(f);
fprintf('written %s\n', fullfile(outDir, 'navic_compare.png'));
end

% ---------------------------------------------------------------- time histories
function figTruth(name, srcPng, dstPng, simDir, runDir, outDir)
setenv('VT_GNSS', 'NavIC');  setenv('VT_RXMODE', '4');  setenv('VT_DATADIR', simDir);
setenv('VT_DATAFILE',       fullfile(simDir, 'IFdata', ['gps_0IF_Fs5MHz_8bit_IQ_Fighter_' name '.bin']));
setenv('VT_DATAFILE_NAVIC', fullfile(simDir, 'IFdata', ['navic_0IF_Fs5MHz_8bit_IQ_Fighter_' name '.bin']));
setenv('VT_IMUFILE',        fullfile(simDir, 'FighterProfiles', ['Fighter_' name '_IMU_200Hz.mat']));
setenv('VT_TRUTHFILE',      fullfile(simDir, 'FighterProfiles', ['Fighter_' name '.csv']));
[file, ~, ~, ~, ~, cmn] = initParameters();
if file.fid > 0, fclose(file.fid); end
if isfield(file, 'fidN') && file.fidN > 0, fclose(file.fidN); end
nav = cell(1,3);  tag = {'CT','VT','UT'};
for m = 1:3
    fn = fullfile(runDir, ['navSol' tag{m} '_' file.fileName '.mat']);
    if exist(fn, 'file'), S = load(fn, ['navSolutions' tag{m}]); nav{m} = S.(['navSolutions' tag{m}]); end
end
S = load(cmn.imuFile, 'imu');
tmp = fullfile(tempdir, ['reportfig_' name]);  if ~exist(tmp, 'dir'), mkdir(tmp); end
compareTruth(nav{1}, nav{2}, nav{3}, file, cmn, tmp, S.imu);
copyfile(fullfile(tmp, srcPng), fullfile(outDir, dstPng));
fprintf('written %s\n', fullfile(outDir, dstPng));
end
