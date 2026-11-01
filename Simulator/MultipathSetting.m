function [mpEnable, mpNumRays] = MultipathSetting(defaultEnable, defaultNumRays)
% -------------------------------------------------------------------------
%                   SoftSim: GPS IF signal simulator
% -------------------------------------------------------------------------
% Prompts for whether to simulate multipath and, if so, how many multipath
% rays to generate per satellite. The statistical multipath model used by
% signalGenerating.m follows:
%   Wang, Y.; Chen, X.; Liu, P. "Statistical Multipath Model Based on
%   Experimental GNSS Data in Static Urban Canyon Environment."
%   Sensors 2018, 18, 1149. https://doi.org/10.3390/s18041149
%
% Optional args (settings.defaultMultipathEnable / .defaultNumMultipathRays)
% skip the interactive prompts for scripted/batch runs. Leave them empty
% ([]) to always be prompted.
%
%   [mpEnable, mpNumRays] = MultipathSetting(defaultEnable, defaultNumRays)

if nargin > 0 && ~isempty(defaultEnable)
    mpEnable = defaultEnable;
else
    fprintf('\nMultipath simulation (statistical model, Wang et al., Sensors 2018, 18, 1149):\n');
    mpEnable = input('Enable multipath? Enter 1 to enable, 0 to disable: ');
    while ~(isscalar(mpEnable) && ismember(mpEnable, [0 1]))
        mpEnable = input('Please enter 1 (enable) or 0 (disable): ');
    end
end

if mpEnable
    if nargin > 1 && ~isempty(defaultNumRays)
        mpNumRays = defaultNumRays;
    else
        mpNumRays = input(['Enter number of multipath rays per satellite ' ...
            '(1-10; the paper found 5 captures most multipath power): ']);
        while ~(isscalar(mpNumRays) && mpNumRays == round(mpNumRays) ...
                && mpNumRays >= 1 && mpNumRays <= 10)
            mpNumRays = input('Please enter an integer between 1 and 10: ');
        end
    end
else
    mpNumRays = 0;
end
