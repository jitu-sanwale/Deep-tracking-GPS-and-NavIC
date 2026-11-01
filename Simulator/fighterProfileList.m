function names = fighterProfileList()
%FIGHTERPROFILELIST Manoeuvre profiles known to generateFighterProfile.m.
%   One per row of the aircraft high-dynamics table, plus a straight-and-
%   level baseline (also the ephemeris reference for the airborne profiles,
%   which all share its initial state).
% Predefined set: a level baseline and seven high-performance manoeuvres.
names = {'baseline_cruise', 'turn', 'loop_pushover', 'break_gonset', 'jinking', ...
         'rolls', 'splitS_immelmann', 'departure_spin'};
end
