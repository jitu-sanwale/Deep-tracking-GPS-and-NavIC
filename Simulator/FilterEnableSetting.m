function filterEn = FilterEnableSetting(defaultEnable)
% -------------------------------------------------------------------------
%                   SoftSim: GPS IF signal simulator
% -------------------------------------------------------------------------
% Prompts for whether to apply the front-end (FE) filter to the generated
% IF signal. Receiver thermal noise (WGN) is added either way - this only
% controls the FIR band-limiting step in signalGenerating.m.
%
% Optional arg (settings.defaultFilterEnable) skips the interactive prompt
% for scripted/batch runs. Leave empty ([]) to always be prompted.
%
%   filterEn = FilterEnableSetting(defaultEnable)

if nargin > 0 && ~isempty(defaultEnable)
    filterEn = defaultEnable;
else
    filterEn = input('\nEnable front-end filter? Enter 1 to enable, 0 to disable: ');
    while ~(isscalar(filterEn) && ismember(filterEn, [0 1]))
        filterEn = input('Please enter 1 (enable) or 0 (disable): ');
    end
end

end
