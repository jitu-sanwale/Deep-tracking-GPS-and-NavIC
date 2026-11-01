function CNoValues = CNoSetting(satList, defaultCNoDB)
% -------------------------------------------------------------------------
%                   SoftSim: GPS IF signal simulator
% Author:
%        Yafeng Li
%    @ Beijing Information Science and Technology University(BISTU)
%    2021. 02. 18
% -------------------------------------------------------------------------
% Optional 2nd arg (settings.defaultCNoDB) skips the interactive prompt,
% so the simulator can run unattended/scripted.
if nargin > 1 && ~isempty(defaultCNoDB)
    if isscalar(defaultCNoDB)
        CNoValues = defaultCNoDB * ones(1,length(satList));
    else
        CNoValues = defaultCNoDB;
    end
    return
end

fprintf('Please set C/No values for all visible satellites in unite of dB-Hz:\n')
fprintf('The visible satellite list is:\n        [');
fprintf('%d ',satList);
disp(']');
fprintf('Format: [C/No1, C/No2, C/No3, ... ] for all %d satellites or one C/No for all satellites. \n',length(satList));
fprintf('Please enter C/No values or press "Ctrl + C" to exit: ');
CNoValues = input('');

while (length(CNoValues) ~= 1) && (length(CNoValues) ~= length(satList))
    fprintf('The number of input must be %d or 1. ',length(satList))
    CNoValues = input('Please enter again: ');
end

if length(CNoValues) == 1
CNoValues = CNoValues * ones(1,length(satList));
end