function S = summarizeFighterRuns(profiles, suffix, runDir, dataDir)
%SUMMARIZEFIGHTERRUNS Cross-profile comparison of the three receiver modes.
%   S = summarizeFighterRuns(profiles, suffix, runDir, dataDir)
%   (defaults: fighterProfileList(), '', Receiver\output, ..\Simulator)
%   suffix '_GPS_NavIC' summarizes the GPS L1 + NavIC L5 runs (results folders
%   ..._Fighter_<name>_GPS_NavIC, output results_FighterSummary_GPS_NavIC) and
%   adds the GPS-NavIC inter-system bias statistics and the NavIC channels used.
%
%   Reads, for each fighter manoeuvre profile,
%     results_gps_0IF_Fs5MHz_8bit_IQ_Fighter_<name>\truth_comparison.mat
%     ..\GPS_Simulator\FighterProfiles\Fighter_<name>_IMU_200Hz.mat (profile stats)
%   and writes to results_FighterSummary\:
%     fighter_summary.csv / .md  - one row per profile x mode: manoeuvre-window
%                                  and whole-run 3D position / velocity error,
%                                  mean channels used; profile dynamics
%     fighter_summary.png        - manoeuvre-window error by profile and mode
%     fighter_summary.json       - the same numbers for a report page

codeDir = fileparts(mfilename('fullpath'));
addpath(fullfile(codeDir, 'common'));
if nargin < 3 || isempty(runDir),  runDir  = fullfile(codeDir, 'output'); end
if nargin < 4 || isempty(dataDir), dataDir = fullfile(codeDir, '..', 'Simulator'); end
if nargin < 1 || isempty(profiles), profiles = fighterProfileList(); end
if nargin < 2, suffix = ''; end
outDir = fullfile(runDir, ['results_FighterSummary' suffix]);
if ~exist(outDir,'dir'), mkdir(outDir); end

modes = {'Classical (DLL/PLL+EKF)', 'Vector (VDCaP+EKF)', 'Ultratight (VDCaP+INS)'};
short = {'CT','VT','UT'};
rows = {};  J = struct('profiles', {{}});
for k = 1:numel(profiles)
    name = profiles{k};
    rf = fullfile(runDir, ['results_gps_0IF_Fs5MHz_8bit_IQ_Fighter_' name suffix], 'truth_comparison.mat');
    pf = fullfile(dataDir, 'FighterProfiles', ['Fighter_' name '_IMU_200Hz.mat']);
    if ~exist(rf,'file') || ~exist(pf,'file')
        fprintf('  %s: no results yet\n', name);  continue;
    end
    R = load(rf);  I = load(pf);  st = I.imu.profile.stats;
    jp = struct('name', name, 'description', I.imu.profile.description, 'model', st.model, ...
        'maxG', st.maxG, 'nRange', [st.minNz st.maxNz], 'onset', st.maxOnsetGps, 'jerk', st.maxJerkGps, ...
        'rollRate', st.maxRollRate, 'alpha', st.alphaRange, 'beta', st.maxBeta, 'hRange', st.hRange, ...
        'VRange', st.VRange, 'inverted', st.invertedFrac, 'window', R.evWin, 'nCh', R.nCh, 'modes', struct());
    for m = 1:numel(modes)
        sm = R.summary(strcmp(R.summary.mode, modes{m}), :);
        if isempty(sm), continue; end
        ev = sm(strcmp(sm.window,'event'),:);  al = sm(strcmp(sm.window,'all'),:);
        isbM = NaN; isbE = NaN; nNav = NaN;
        if isfield(R,'isbSummary') && ~isempty(R.isbSummary)
            ib = R.isbSummary(strcmp(R.isbSummary.mode, modes{m}) & strcmp(R.isbSummary.window,'event'), :);
            if ~isempty(ib), isbM = ib.isb_mean_m; isbE = ib.isb_err_rms_m; nNav = ib.navicUsed_mean; end
        end
        rows(end+1,:) = {name, short{m}, ev.pos3D_rms_m, ev.pos3D_max_m, ev.vel3D_rms_mps, ev.vel3D_max_mps, ...
            ev.chUsed_mean, al.pos3D_rms_m, al.vel3D_rms_mps, st.maxG, st.maxOnsetGps, st.maxRollRate, ...
            100*st.invertedFrac, isbM, isbE, nNav}; %#ok<AGROW>
        jp.modes.(short{m}) = struct('evPosRms', ev.pos3D_rms_m, 'evPosMax', ev.pos3D_max_m, ...
            'evVelRms', ev.vel3D_rms_mps, 'evVelMax', ev.vel3D_max_mps, 'evCh', ev.chUsed_mean, ...
            'allPosRms', al.pos3D_rms_m, 'allVelRms', al.vel3D_rms_mps);
        ts = R.errTS(strcmp({R.errTS.name}, modes{m}));
        if ~isempty(ts)
            i = 1:max(1, round(numel(ts.t)/600)):numel(ts.t);     % <= ~600 points per series
            jp.modes.(short{m}).ts = struct('t', round(ts.t(i)',2), 'e3', round(vecnorm(ts.eP(i,:),2,2)',3), ...
                'v3', round(vecnorm(ts.eV(i,:),2,2)',4), 'used', ts.used(i)');
        end
    end
    J.profiles{end+1} = jp;
end
if isempty(rows), error('summarizeFighterRuns: no results found.'); end
S = cell2table(rows, 'VariableNames', {'profile','mode','man_pos_rms_m','man_pos_max_m','man_vel_rms_mps', ...
    'man_vel_max_mps','man_chUsed','all_pos_rms_m','all_vel_rms_mps','maxG','onset_gps','rollrate_dps','inverted_pct', ...
    'man_isb_mean_m','man_isb_err_rms_m','man_navicUsed'});
writetable(S, fullfile(outDir, 'fighter_summary.csv'));
fid = fopen(fullfile(outDir, 'fighter_summary.json'), 'w');  fprintf(fid, '%s', jsonencode(J));  fclose(fid);

% markdown table
fid = fopen(fullfile(outDir, 'fighter_summary.md'), 'w');
fprintf(fid, '| Profile | Mode | Manoeuvre pos rms / max (m) | Manoeuvre vel rms / max (m/s) | Channels used | Whole-run pos / vel rms |\n');
fprintf(fid, '|---|---|---|---|---|---|\n');
fm = @(x, f) fmtDiv(x, f);
for i = 1:height(S)
    fprintf(fid, '| %s | %s | %s / %s | %s / %s | %.1f | %s m / %s m/s |\n', strrep(S.profile{i},'_','\_'), S.mode{i}, ...
        fm(S.man_pos_rms_m(i),'%.1f'), fm(S.man_pos_max_m(i),'%.1f'), fm(S.man_vel_rms_mps(i),'%.2f'), ...
        fm(S.man_vel_max_mps(i),'%.2f'), S.man_chUsed(i), fm(S.all_pos_rms_m(i),'%.1f'), fm(S.all_vel_rms_mps(i),'%.2f'));
end
fclose(fid);

% figure: manoeuvre-window rms error, grouped by profile (log scale)
pr = unique(S.profile, 'stable');
col = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478];
f = figure('Color','w','Position',[60 60 1300 760],'Visible','off');
tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact');
for q = 1:2
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
    set(ax, 'YScale','log', 'GridAlpha',0.12, 'FontSize',9, 'TickDir','out');
    M = nan(numel(pr), 3);
    for i = 1:numel(pr)
        for m = 1:3
            r = strcmp(S.profile, pr{i}) & strcmp(S.mode, short{m});
            if any(r)
                if q == 1, M(i,m) = S.man_pos_rms_m(r); else, M(i,m) = S.man_vel_rms_mps(r); end
            end
        end
    end
    b = bar(ax, max(M, 1e-3), 'grouped', 'EdgeColor', 'none');
    for m = 1:3, b(m).FaceColor = col(m,:); end
    set(ax, 'XTick', 1:numel(pr), 'XTickLabel', strrep(pr,'_','\_'));
    if q == 1
        ylabel(ax, '3D position rms (m)');
        title(ax, 'Error vs truth during each manoeuvre window (log scale)');
        legend(ax, modes, 'Location', 'northwest', 'Box', 'off');
    else
        ylabel(ax, '3D velocity rms (m/s)');
    end
end
exportgraphics(f, fullfile(outDir, 'fighter_summary.png'), 'Resolution', 150);  close(f);
disp(S);
fprintf('Fighter summary written to %s\n', outDir);
end

function s = fmtDiv(x, f)
% a numerically diverged solution (non-finite or > 1e6) is reported as such
if ~isfinite(x) || abs(x) > 1e6, s = 'diverged'; else, s = sprintf(f, x); end
end
