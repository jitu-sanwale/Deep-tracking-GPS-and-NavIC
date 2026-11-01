function modes = rxModeSetting(defaultMode)
%Purpose:
%   Select the receiver mode(s) to run:
%     1 - classical:  scalar DLL/PLL per channel + 8-state EKF   (trackingCT_POS)
%     2 - vector:     VDCaP (VDLL + Doppler-aided PLL) + 12-state
%                     Singer PVA EKF, no INS        (trackingUT_POS 'none')
%     3 - ultratight: VDCaP (VDLL + Doppler-aided PLL) + 18-state
%                     EKF + 6-DOF INS                             (trackingUT_POS)
%     4 - all three (default)
%   Modes 2 and 3 are handed over from the classical solution, so classical
%   tracking always runs; it is reported only when selected.
%   defaultMode (or the environment variable VT_RXMODE) skips the prompt.
%Outputs:
%   modes - logical [CT VT UT]
if nargin < 1, defaultMode = []; end
if isempty(defaultMode) && ~isempty(getenv('VT_RXMODE'))
    defaultMode = str2double(getenv('VT_RXMODE'));
end
if ~isempty(defaultMode)
    choice = defaultMode;
else
    fprintf('\nReceiver mode:\n');
    fprintf('  1) Classical  - DLL/PLL + EKF\n');
    fprintf('  2) Vector     - VDLL + EKF\n');
    fprintf('  3) Ultratight - VDCaP + 17-state EKF + 6-DOF INS (iNAT-FSLG-01)\n');
    fprintf('  4) All three, compared against the truth trajectory (default)\n');
    try
        choice = input('Select [4]: ');
    catch
        fprintf('  (no interactive input available - using the default)\n');
        choice = 4;
    end
    if isempty(choice), choice = 4; end
end
switch choice
    case 1, modes = [true false false];
    case 2, modes = [false true false];
    case 3, modes = [false false true];
    case 4, modes = [true true true];
    otherwise, error('rxModeSetting: invalid mode %g (expected 1-4).', choice);
end
names = {'classical','vector','ultratight'};
fprintf('Receiver mode(s): %s\n', strjoin(names(modes), ', '));
end
