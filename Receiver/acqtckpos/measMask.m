function [rows, used] = measMask(cn0, el, prInno, msIndex, track)
%Purpose:
%   Select the channels used in the EKF update (excluded channels keep being
%   tracked; VT keeps steering them from the navigation solution):
%     - C/N0 below track.cn0Mask (lost or very weak signal)
%     - elevation below track.elMask
%     - |pseudorange innovation| above track.prGate (after track.gateStart ms),
%       applied only while at least 4 channels would remain
%   rows - measurement rows to keep ([pseudorange rows, pseudorange-rate rows])
n    = numel(cn0);
used = true(1,n);
used(cn0 < track.cn0Mask) = false;
if ~isempty(el), used(el(1:n) < track.elMask) = false; end
if msIndex > track.gateStart && ~isempty(prInno)
    g = used & (abs(prInno(:)') <= track.prGate);
    if sum(g) >= 4, used = g; end
end
if sum(used) < 4          % never starve the filter
    [~, ord] = sort(cn0, 'descend');
    used(ord(1:min(4,n))) = true;
end
rows = [find(used), find(used)+n];
end
