function summary = compareTruth(navCT, navVT, navUT, file, cmn, outDir, imu)
%Purpose:
%   Compare the classical (CT), vector (VT) and ultratight (UT) navigation
%   solutions against the truth trajectory the IF file was generated from.
%   Truth = ECEF cubic spline of cmn.truthFile (trajectorySplineECEF.m in
%   Receiver\common), the same spline that produced the simulator's 1 kHz input,
%   with exact spline-derivative velocity. Epoch times are data-file times
%   (sample index / Fs), which equal trajectory time in the simulation.
%   Any of navCT/navVT/navUT may be empty.
%Outputs:
%   summary - table of error statistics (also written to outDir)
if ~exist(outDir,'dir'), mkdir(outDir); end
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'common'));
traj = trajectorySplineECEF(cmn.truthFile);

T = 1e-3;
sol = struct('name',{},'t',{},'pos',{},'vel',{},'used',{});
usedOf = @(nav) sum(nav.usedCh, 2);
if ~isempty(navCT)
    sol(end+1) = struct('name','Classical (DLL/PLL+EKF)', 't',(file.skip + (1:size(navCT.usrPos,1))')*T, ...
                        'pos',navCT.usrPos, 'vel',navCT.usrVel, 'used',usedOf(navCT));
end
if ~isempty(navVT)
    sol(end+1) = struct('name','Vector (VDCaP+EKF)', 't',navVT.tData(:), 'pos',navVT.usrPos, 'vel',navVT.usrVel, ...
                        'used',usedOf(navVT));
end
if ~isempty(navUT)
    sol(end+1) = struct('name','Ultratight (VDCaP+INS)', 't',navUT.tData(:), 'pos',navUT.usrPos, 'vel',navUT.usrVel, ...
                        'used',usedOf(navUT));
end

% Reference palette slots 1-3 (blue, orange, aqua) + line style as a second channel
col = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478];
sty = {'-','--','-'};
lbl = {'East','North','Up'};
tEv = [25 27];                                % spliced jerk event (ProjectileTrajectory_JerkSplice)
evName = '120 g/s jerk event';
if nargin >= 7 && ~isempty(imu) && isfield(imu,'profile')
    tEv = imu.profile.maneuverWin;  evName = ['manoeuvre (' strrep(imu.profile.name,'_','\_') ')'];
end

%% errors in local ENU at the truth position
for k = 1:numel(sol)
    t  = sol(k).t;  ok = t >= traj.t0 & t <= traj.t1;
    t  = t(ok);
    rT = ppval(traj.ppPos, t')';  vT = ppval(traj.ppVel, t')';
    dr = sol(k).pos(ok,:) - rT;   dv = sol(k).vel(ok,:) - vT;
    eP = zeros(numel(t),3); eV = zeros(numel(t),3);
    llh = zeros(numel(t),2);
    for i = 1:numel(t)
        q = xyz2llh(rT(i,:)); llh(i,:) = q(1:2);
    end
    sL = sin(llh(:,1)); cL = cos(llh(:,1)); sl = sin(llh(:,2)); cl = cos(llh(:,2));
    rot = @(d) [-sl.*d(:,1) + cl.*d(:,2), ...
                -sL.*cl.*d(:,1) - sL.*sl.*d(:,2) + cL.*d(:,3), ...
                 cL.*cl.*d(:,1) + cL.*sl.*d(:,2) + sL.*d(:,3)];
    eP = rot(dr);  eV = rot(dv);
    sol(k).t = t; sol(k).eP = eP; sol(k).eV = eV;  sol(k).used = sol(k).used(ok);  sol(k).okMask = ok;
    sol(k).e3 = vecnorm(eP,2,2); sol(k).v3 = vecnorm(eV,2,2);
end

%% statistics: whole common window, jerk event, after the event
tStart = file.skip*T + file.skiptimeVT*T + 1;     % 1 s after the VT/UT handover
tEnd   = min(arrayfun(@(s) s.t(end), sol));
% every window is clipped to the interval common to all receivers (after the
% VT/UT handover and before the earliest end of record), so CT, VT and UT are
% compared over the same epochs
win = {'all', [tStart tEnd]; 'event', [max(tEv(1), tStart) min(tEv(2)+1, tEnd)]; 'post', [tEv(2)+1 tEnd]};
rows = {};
mx = @(x) max([x(:); NaN], [], 'omitnan');            % NaN for an empty window
rm = @(x) sqrt(mean(x(:).^2, 'omitnan'));             % a diverged (NaN) solution is skipped
for k = 1:numel(sol)
    for w = 1:size(win,1)
        m = sol(k).t >= win{w,2}(1) & sol(k).t <= win{w,2}(2);
        e = sol(k).eP(m,:);
        rows(end+1,:) = {sol(k).name, win{w,1}, rm(sol(k).e3(m)), mx(sol(k).e3(m)), ...
            rm(vecnorm(e(:,1:2),2,2)), rm(e(:,3)), rm(sol(k).v3(m)), mx(sol(k).v3(m)), mean(sol(k).used(m))}; %#ok<AGROW>
    end
end
summary = cell2table(rows, 'VariableNames', {'mode','window','pos3D_rms_m','pos3D_max_m', ...
    'horiz_rms_m','up_rms_m','vel3D_rms_mps','vel3D_max_mps','chUsed_mean'});
writetable(summary, fullfile(outDir, 'truth_comparison_summary.csv'));
errTS = struct('name',{},'t',{},'eP',{},'eV',{},'used',{});   % 10 Hz series for cross-run summaries
for k = 1:numel(sol)
    i = 1:max(1, round(0.1/median(diff(sol(k).t)))):numel(sol(k).t);
    errTS(k) = struct('name', sol(k).name, 't', sol(k).t(i), 'eP', sol(k).eP(i,:), 'eV', sol(k).eV(i,:), 'used', sol(k).used(i));
end
evWin = tEv; nCh = 0;
if ~isempty(navCT), nCh = size(navCT.usedCh,2); end

% GPS + NavIC: estimated GPS-NavIC inter-system bias against its truth (the
% receiver L5-L1 delay minus the error of the broadcast IRNSS-GPS offset,
% written by the simulator next to the NavIC IF file) and NavIC channels used
isbSummary = table();
navs = {navCT, navVT, navUT};  names = {'Classical (DLL/PLL+EKF)','Vector (VDCaP+EKF)','Ultratight (VDCaP+INS)'};
tf = '';
if isfield(file,'fileRouteN'), [d_, b_] = fileparts(file.fileRouteN); tf = fullfile(d_, [b_ '_time.mat']); end
if ~isempty(tf) && exist(tf,'file') && ~isempty(navCT) && isfield(navCT,'sv') && any(navCT.sv > 100)
    tt = load(tf);  isN = navCT.sv > 100;
    rowsI = {};
    for k = 1:3
        nv = navs{k};
        if isempty(nv) || ~isfield(nv,'isb'), continue; end
        s_ = sol(strcmp({sol.name}, names{k}));
        tI = s_.t;  isb = nv.isb(:);  isb = isb(s_.okMask);             % rows kept by the truth window
        uN = sum(nv.usedCh(s_.okMask, isN), 2);
        for w = 1:size(win,1)
            m = tI >= win{w,2}(1) & tI <= win{w,2}(2);
            rowsI(end+1,:) = {names{k}, win{w,1}, tt.isbTruth_m, mean(isb(m),'omitnan'), ...
                rm(isb(m) - tt.isbTruth_m), mean(uN(m))}; %#ok<AGROW>
        end
    end
    isbSummary = cell2table(rowsI, 'VariableNames', {'mode','window','isb_truth_m','isb_mean_m','isb_err_rms_m','navicUsed_mean'});
    writetable(isbSummary, fullfile(outDir, 'isb_summary.csv'));
    fprintf('\nGPS-NavIC inter-system bias (truth %.3f m = c x (L5 delay %.1f ns - offset error %.1f ns)):\n', ...
        tt.isbTruth_m, tt.rxL5Delay*1e9, tt.dtErr*1e9);
    disp(isbSummary);
end
save(fullfile(outDir, 'truth_comparison.mat'), 'summary', 'errTS', 'evWin', 'nCh', 'win', 'isbSummary');
fprintf('\nError vs truth (windows: all = %.1f-%.1f s, event = %g-%g s, post = %g-%.1f s)\n', ...
    tStart, tEnd, win{2,2}(1), win{2,2}(2), tEv(2)+1, tEnd);
disp(summary);

% 10 Hz error time series per mode
for k = 1:numel(sol)
    step = max(1, round(0.1/median(diff(sol(k).t))));
    i = 1:step:numel(sol(k).t);
    tag = lower(strtok(sol(k).name));
    writematrix([sol(k).t(i) sol(k).eP(i,:) sol(k).eV(i,:)], fullfile(outDir, ['err_' tag '_10Hz.csv']));
end

%% figures
f = newFig('Position error'); tl = tiledlayout(f,3,1,'TileSpacing','compact');
for a = 1:3
    ax = nexttile(tl); styleAx(ax); shadeEvent(ax, tEv);
    for k = 1:numel(sol), plot(ax, sol(k).t, sol(k).eP(:,a), sty{k}, 'Color', col(k,:), 'LineWidth', 2.2); end
    ylabel(ax, [lbl{a} ' error (m)']);
    if a == 1, title(ax, ['Position error vs truth (ENU) - shaded: ' evName]);
        legend(ax, {sol.name}, 'Location','northwest', 'Box','on', 'TextColor','k', 'FontSize',14); end
end
xlabel(ax, 'Time (s)'); saveFig(f, fullfile(outDir, 'T1_position_error_ENU.png'));

f = newFig('Velocity error'); tl = tiledlayout(f,3,1,'TileSpacing','compact');
for a = 1:3
    ax = nexttile(tl); styleAx(ax); shadeEvent(ax, tEv);
    for k = 1:numel(sol), plot(ax, sol(k).t, sol(k).eV(:,a), sty{k}, 'Color', col(k,:), 'LineWidth', 2.2); end
    ylabel(ax, ['V_' lbl{a}(1) ' error (m/s)']);
    if a == 1, title(ax, 'Velocity error vs truth (ENU)');
        legend(ax, {sol.name}, 'Location','northwest', 'Box','on', 'TextColor','k', 'FontSize',14); end
end
xlabel(ax, 'Time (s)'); saveFig(f, fullfile(outDir, 'T2_velocity_error_ENU.png'));

f = newFig('3D error'); tl = tiledlayout(f,2,1,'TileSpacing','compact');
ax = nexttile(tl); styleAx(ax); shadeEvent(ax, tEv); set(ax,'YScale','log');
for k = 1:numel(sol), plot(ax, sol(k).t, max(sol(k).e3,1e-3), sty{k}, 'Color', col(k,:), 'LineWidth', 2.2); end
ylabel(ax, '3D position error (m)'); title(ax, '3D error magnitude (log scale)');
legend(ax, {sol.name}, 'Location','northwest', 'Box','on', 'TextColor','k', 'FontSize',14);
ax = nexttile(tl); styleAx(ax); shadeEvent(ax, tEv); set(ax,'YScale','log');
for k = 1:numel(sol), plot(ax, sol(k).t, max(sol(k).v3,1e-4), sty{k}, 'Color', col(k,:), 'LineWidth', 2.2); end
ylabel(ax, '3D velocity error (m/s)'); xlabel(ax, 'Time (s)');
saveFig(f, fullfile(outDir, 'T3_error_magnitude.png'));

% ultratight diagnostics
if ~isempty(navUT)
    t = navUT.tData;
    if nargin >= 7 && ~isempty(imu)
        % small-angle attitude error from the DCMs (robust through vertical /
        % inverted flight, where Euler-angle differences jump by 180 deg)
        attT = interp1(imu.truth.t, unwrap(imu.truth.attRad), t, 'nearest', 'extrap');
        dAtt = zeros(numel(t),3);
        for i = 1:numel(t)
            dC = cbnE(navUT.attEuler(i,:)*pi/180) * cbnE(attT(i,:))';
            dAtt(i,:) = [dC(3,2) dC(1,3) dC(2,1)]*180/pi;
        end
        attLbl = {'about N','about E','about D'};
    else
        attT = ppval(traj.ppAtt, t')'*180/pi;
        dAtt = navUT.attEuler - attT;  dAtt(:,3) = mod(dAtt(:,3)+180,360)-180;
        attLbl = {'Roll','Pitch','Yaw'};
    end
    f = newFig('UT diagnostics'); tl = tiledlayout(f,4,1,'TileSpacing','compact');
    ax = nexttile(tl); styleAx(ax); shadeEvent(ax, tEv);
    plot(ax, t, dAtt(:,1:2), 'LineWidth', 1.8); plot(ax, t, dAtt(:,3), 'LineWidth', 2.2);
    colororder(ax, col); ylabel(ax, 'Attitude error (deg)');
    legend(ax, attLbl, 'Location','northeast', 'Box','on', 'TextColor','k', 'FontSize',14);
    title(ax, 'Ultratight: attitude, IMU bias and CaP loop status');
    if nargin >= 7 && ~isempty(imu)
        baT = interp1(imu.t, imu.biasA, t)/9.80665e-6;
        ax = nexttile(tl); styleAx(ax); shadeEvent(ax, tEv); colororder(ax, col);
        plot(ax, t, navUT.ba/9.80665e-6, 'LineWidth', 2.2); plot(ax, t, baT, ':', 'LineWidth', 1.8);
        ylabel(ax, 'Accel bias (\mug)'); legend(ax, {'x est','y est','z est','x true','y true','z true'}, ...
            'Location','eastoutside', 'Box','on', 'TextColor','k', 'FontSize',14);
    end
    ax = nexttile(tl); styleAx(ax); shadeEvent(ax, tEv);
    plot(ax, t, min(navUT.PLI,[],2), 'Color', col(1,:), 'LineWidth', 1.8);
    plot(ax, t, mean(navUT.PLI,2), 'Color', col(2,:), 'LineWidth', 2.2);
    yline(ax, 0.9, ':', 'reset threshold', 'LineWidth', 1.6, 'FontSize', 13, 'LabelHorizontalAlignment','left'); ylim(ax, [-0.2 1.05]);
    ylabel(ax, 'PLI'); legend(ax, {'min over channels','mean'}, 'Location','southwest', 'Box','on', 'TextColor','k', 'FontSize',14);
    ax = nexttile(tl); styleAx(ax); shadeEvent(ax, tEv);
    stairs(ax, t, sum(navUT.usedCh,2), 'Color', col(1,:), 'LineWidth', 2.2);
    stairs(ax, t, sum(navUT.reset,2), 'Color', col(2,:), 'LineWidth', 1.8);
    ylabel(ax, 'Channels'); xlabel(ax, 'Time (s)');
    legend(ax, {'used in EKF','CaP resets'}, 'Location','west', 'Box','on', 'TextColor','k', 'FontSize',14);
    saveFig(f, fullfile(outDir, 'T4_ultratight_diagnostics.png'));
end
fprintf('Truth comparison saved in %s\n', outDir);
end

function C = cbnE(e)
% body -> NED for Euler [roll pitch yaw] (rad)
sp = sin(e(1)); cp = cos(e(1)); st = sin(e(2)); ct = cos(e(2)); ss = sin(e(3)); cs = cos(e(3));
C = [ct*cs, ct*ss, -st; -cp*ss+sp*st*cs, cp*cs+sp*st*ss, sp*ct; sp*ss+cp*st*cs, -sp*cs+cp*st*ss, cp*ct]';
end

function f = newFig(name)
f = figure('Name', name, 'Color', 'w', 'Position', [60 60 820 760], 'Visible', 'off');
end
function styleAx(ax)
hold(ax,'on'); grid(ax,'on'); box(ax,'on');
set(ax, 'GridAlpha', 0.15, 'GridColor', [0.3 0.3 0.3], 'XColor', 'k', ...
    'YColor', 'k', 'FontSize', 15, 'TickDir', 'out', 'LineWidth', 1.4);
ax.Title.Color = 'k';  ax.Title.FontWeight = 'normal';
end
function shadeEvent(ax, tEv)
xregion(ax, tEv(1), tEv(2), 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.12, 'HandleVisibility', 'off');
end
function saveFig(f, fn)
exportgraphics(f, fn, 'Resolution', 150);
close(f);
end
