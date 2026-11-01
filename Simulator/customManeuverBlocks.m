function B = customManeuverBlocks()
%CUSTOMMANEUVERBLOCKS Building blocks for a user-defined fighter manoeuvre.
%   B = customManeuverBlocks() returns a struct array with fields
%     name     - block keyword used in opt.customPhases ({'turn', 7, 10, 3, 1}, ...)
%     params   - parameter names, in order
%     defaults - default values (used when a parameter is omitted)
%     lo, hi   - accepted range of each parameter
%     help     - one-line description
%   The blocks are flown in sequence by the nonlinear 6-DOF F-16 under the
%   dynamic-inversion autopilot (fighterPhases.m, case 'custom'). Load-factor
%   commands are capped by the available lift, as in the predefined profiles, so a
%   demand the aircraft cannot meet shows up as a lower g, not a stall.

B = struct('name', {}, 'params', {}, 'defaults', {}, 'lo', {}, 'hi', {}, 'help', {});
B(end+1) = blk('cruise', {'duration_s'},                              10,           1,            60,          'straight and level flight');
B(end+1) = blk('turn',   {'peak_g','hold_s','onset_gps','dir'}, [7 10 3 1], [1.5 0 0.5 -1], [9 30 15 1], 'level turn (dir +1 right, -1 left): g onset, hold, release, roll-out');
B(end+1) = blk('bank',   {'bank_deg','rate_dps','hold_s'},             [60 180 3],   [-90 10 0],   [90 300 30], 'roll to a bank angle and hold it (level flight)');
B(end+1) = blk('roll',   {'angle_deg','rate_dps'},                     [360 270],    [-1080 30],   [1080 300],  'aileron roll through an angle at n = 1 (+ right, - left)');
B(end+1) = blk('pull',   {'n_g','duration_s'},                         [5 2],        [-3 0.2],     [9 20],      'wings-level pull (n > 1) or push (n < 1) for a duration');
B(end+1) = blk('loop',   {'peak_g'},                                   6,            3,            9,           'full loop at up to peak_g, then 4 s recovery');
B(end+1) = blk('unload', {'n_g','duration_s'},                         [1.5 0.8],    [-1 0.2],     [3 5],       'hold the bank and unload to n_g (before a roll-out)');
end

function b = blk(name, params, defaults, lo, hi, help)
b = struct('name', name, 'params', {params}, 'defaults', defaults, 'lo', lo, 'hi', hi, 'help', help);
end
