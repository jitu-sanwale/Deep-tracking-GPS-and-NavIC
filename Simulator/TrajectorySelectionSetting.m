function [trajFile, singleSatEnable] = TrajectorySelectionSetting(defaultTrajFile, defaultSingleSatEnable)
% -------------------------------------------------------------------------
%                   SoftSim: GPS IF signal simulator
% -------------------------------------------------------------------------
% Prompts for which ground-truth trajectory CSV to load (see
% signalGenerating.m's "Load trajectory" section for the expected column
% format: Time [sec], Latitude [deg], Longitude [deg], Height [meter],
% Roll Angle [deg], Pitch Angle [deg], Heading [deg]) and whether to
% restrict signal generation to a single satellite.
%
% The trajectory menu lists every *.csv found in the current folder, plus
% two options to (re)generate a stress-test trajectory:
%
%   - High-dynamics jerk-doublet trajectory (see
%     generateHighDynamicsTrajectory.m): jerk values [0 20 40 80 100 120]
%     m/s^3, each a doublet maneuver (jerk +J then -J, Th=1s each)
%     starting at its own timestamp along one continuous, single-file
%     trajectory, with cruise segments in between so the loop can
%     re-settle before the next event. Choosing this option (re)writes
%     HighDynamicsTrajectory.csv in the current folder and selects it.
%
%   - Fig.7-style high-dynamics trajectory (see
%     generateFig7Trajectory.m): a 60 s trajectory containing one 7 s
%     maneuver reproducing the acceleration/jerk shape of Won, Pany &
%     Eissfeller (IEEE TAES 48(4), 2012), Fig. 7, scaled to a
%     user-prompted peak jerk (default 120 g/s, vs. the paper's own
%     100 g/s). The maneuver is placed starting at a user-prompted time
%     (default 20 s, i.e. spanning [20, 27] s), with constant-velocity
%     cruise before and after so the full 60 s run demonstrates both
%     tracking-loop performance through the maneuver and (with
%     single-satellite mode left OFF, so a full multi-satellite position
%     fix is possible) position error against the known ground truth.
%     Choosing this option (re)writes Fig7Trajectory.csv in the current
%     folder and selects it.
%
% Single-satellite mode restricts the visible-satellite list (computed in
% signalGenerating.m via getVisibleSat.m) down to the single
% highest-elevation satellite, right after it is determined - this
% isolates one tracking channel from the confound of other satellites'
% independently-varying geometry/dynamics. It is offered for any
% trajectory; for the Fig.7-style trajectory the default suggestion is
% OFF (0), unlike the jerk-doublet trajectory, since demonstrating
% position error needs the normal multi-satellite visible set.
%
% Optional args (settings.defaultTrajFile / .defaultSingleSatEnable) skip
% the interactive prompts for scripted/batch runs. Leave empty ([]) to
% always be prompted.
%
%   [trajFile, singleSatEnable] = TrajectorySelectionSetting(defaultTrajFile, defaultSingleSatEnable)

HIGH_DYN_FILE  = 'HighDynamicsTrajectory.csv';
HIGH_DYN_JERKS = [0 20 40 80 100 120];   % [m/s^3]

FIG7_FILE            = 'Fig7Trajectory.csv';
FIG7_DEFAULT_START   = 20;    % [s] maneuver start (spans [start, start+7] s)
FIG7_DEFAULT_JERK_GPS = 120;  % [g/s] peak jerk (paper's own scenario: 100 g/s)

if nargin > 0 && ~isempty(defaultTrajFile)
    trajFile = defaultTrajFile;
    % Scripted/batch use of the Fig.7 trajectory: regenerate it with
    % whatever parameters the caller set on the base workspace / settings
    % struct isn't available here, so a batch driver that wants a
    % Fig.7Trajectory.csv should call generateFig7Trajectory(...) itself
    % before invoking signalGenerating with settings.defaultTrajFile set
    % to FIG7_FILE - this function only regenerates it interactively,
    % below.
else
    csvList  = dir('*.csv');
    csvNames = {csvList.name};
    fprintf('\nGround-truth trajectory selection:\n');
    for k = 1:numel(csvNames)
        fprintf('  %d) %s\n', k, csvNames{k});
    end
    genOpt1 = numel(csvNames) + 1;
    genOpt2 = numel(csvNames) + 2;
    fprintf('  %d) Generate high-dynamics jerk-doublet trajectory (%s)\n', genOpt1, HIGH_DYN_FILE);
    fprintf('     Jerk values %s m/s^3 as successive doublet maneuvers along one\n', mat2str(HIGH_DYN_JERKS));
    fprintf('     continuous trajectory, for single-satellite tracking-loop stress testing.\n');
    fprintf('  %d) Generate Fig.7-style high-dynamics trajectory (%s)\n', genOpt2, FIG7_FILE);
    fprintf('     A 60 s trajectory with one 7 s maneuver (default: starts at %g s, peak jerk\n', FIG7_DEFAULT_START);
    fprintf('     %g g/s) reproducing Won/Pany/Eissfeller (2012) Fig. 7''s acceleration/jerk\n', FIG7_DEFAULT_JERK_GPS);
    fprintf('     shape, for demonstrating BOTH tracking-loop performance and position error.\n');
    choice = input(sprintf('Select trajectory [1-%d]: ', genOpt2));
    while ~(isscalar(choice) && choice == round(choice) && choice >= 1 && choice <= genOpt2)
        choice = input(sprintf('Please enter an integer between 1 and %d: ', genOpt2));
    end
    if choice == genOpt1
        generateHighDynamicsTrajectory(HIGH_DYN_FILE, HIGH_DYN_JERKS);
        trajFile = HIGH_DYN_FILE;
    elseif choice == genOpt2
        maneuverStart = input(sprintf('  Maneuver start time [s] (default %g): ', FIG7_DEFAULT_START));
        if isempty(maneuverStart), maneuverStart = FIG7_DEFAULT_START; end
        jerkPeakGps = input(sprintf('  Peak jerk [g/s] (default %g, paper=100): ', FIG7_DEFAULT_JERK_GPS));
        if isempty(jerkPeakGps), jerkPeakGps = FIG7_DEFAULT_JERK_GPS; end
        generateFig7Trajectory(FIG7_FILE, maneuverStart, jerkPeakGps);
        trajFile = FIG7_FILE;
    else
        trajFile = csvNames{choice};
    end
end

isHighDyn = strcmpi(trajFile, HIGH_DYN_FILE);
isFig7    = strcmpi(trajFile, FIG7_FILE);
if nargin > 1 && ~isempty(defaultSingleSatEnable)
    singleSatEnable = defaultSingleSatEnable;
else
    if isHighDyn
        defaultAns = 1;   % single-channel Doppler stress test
    elseif isFig7
        defaultAns = 0;   % needs the full multi-satellite set for position error
    else
        defaultAns = 0;
    end
    singleSatEnable = input(sprintf(...
        'Restrict to a single (highest-elevation) satellite? Enter 1 to enable, 0 to disable [default %d]: ', ...
        defaultAns));
    if isempty(singleSatEnable)
        singleSatEnable = defaultAns;
    end
    while ~(isscalar(singleSatEnable) && ismember(singleSatEnable, [0 1]))
        singleSatEnable = input('Please enter 1 (enable) or 0 (disable): ');
    end
end
