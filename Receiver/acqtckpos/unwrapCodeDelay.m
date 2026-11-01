function cd = unwrapCodeDelay(cdPrev, cdMod, period)
%Purpose:
%   Keep the code delay (in samples) continuous across the 1-ms boundary so
%   that pseudoranges do not jump by 1 ms (~300 km) when the code phase wraps.
d  = cdMod - mod(cdPrev, period);
d  = d - period*round(d/period);
cd = cdPrev + d;
end
