function rangeG = AccelRangeSetting(defaultRangeG)
% -------------------------------------------------------------------------
% Prompts for the accelerometer measurement range (saturation limit) used by
% generateIMU_iNAT_FSLG01.m, same prompt pattern as the other *Setting.m
% files of this simulator.
%
%   1) Very high (default) - +/-10000 g, i.e. effectively no clipping, so the
%      IMU error model is the datasheet one minus the range limit. The
%      ProjectileTrajectory_JerkSplice flight peaks at ~60 g during the
%      spliced jerk event, three times the real sensor's range.
%   2) iNAT-FSLG-01 datasheet - +/-20 g. The accelerometers clip during the
%      jerk event (t ~ 25-27 s), exactly as the real unit would.
%   3) Custom value [g].
%
% Pass defaultRangeG (e.g. settings.defaultAccelRangeG) to skip the prompt
% for scripted runs; leave it empty ([]) to be prompted.
%
%   rangeG = AccelRangeSetting(defaultRangeG)
% -------------------------------------------------------------------------
VERY_HIGH_G  = 1e4;
DATASHEET_G  = 20;

if nargin > 0 && ~isempty(defaultRangeG)
    rangeG = defaultRangeG;
else
    fprintf('\nAccelerometer measurement range (saturation):\n');
    fprintf('  1) Very high, +/-%g g - effectively unlimited (default)\n', VERY_HIGH_G);
    fprintf('  2) iNAT-FSLG-01 datasheet, +/-%g g (clips during the ~60 g jerk event)\n', DATASHEET_G);
    fprintf('  3) Custom value\n');
    try
        choice = input('Select [1]: ');
    catch
        % no interactive input (e.g. matlab -batch): take the default
        fprintf('  (no interactive input available - using the default)\n');
        choice = 1;
    end
    if isempty(choice), choice = 1; end
    switch choice
        case 1, rangeG = VERY_HIGH_G;
        case 2, rangeG = DATASHEET_G;
        case 3
            rangeG = input('Accelerometer range [g]: ');
            if isempty(rangeG), rangeG = VERY_HIGH_G; end
        otherwise
            error('AccelRangeSetting: invalid selection %g (expected 1, 2 or 3).', choice);
    end
end
if ~isscalar(rangeG) || ~(rangeG > 0)
    error('AccelRangeSetting: range must be a positive scalar in g (got %s).', mat2str(rangeG));
end
fprintf('Accelerometer range: +/-%g g.\n', rangeG);
end
