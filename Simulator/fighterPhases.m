function P = fighterPhases(name, opt)
%FIGHTERPHASES Manoeuvre definitions for generateFighterProfile.m.
%   P = fighterPhases(name, opt)
%
%   A profile is a sequence of phases. Each phase has
%     cmd(tau, s, s0) -> [p_ref (rad/s), n_ref (g), dV/dt (m/s^2), onGround,
%                         n_direct (g, point mass only), surface override]
%     done(tau, s, s0) -> true when the phase ends
%   tau is the time since the phase started, s the current state, s0 the
%   state at the phase start (fields V h phi theta gamma psi qint pint nAct
%   alpha w upright). The surface override ([elevator aileron rudder
%   throttle], NaN = autopilot) is used by the 6-DOF model only.
%
%   Airborne profiles start at 250 m/s, 5000 m, heading East, and share
%   their first 12 s of straight and level flight (so they can share the
%   baseline's decoded ephemeris); they fly the nonlinear 6-DOF F-16 (NASA
%   TP-1538 aerodynamics, F16Model.jl port). The ground-contact profiles use the
%   point-mass model.
%   Fighter-level targets: 9 g / 12 g/s onset breaks, 270-300 deg/s rolls,
%   negative-g pushes, inverted flight, departure, 4.5 g catapult, 4-5 g
%   arrested landing and touchdown spikes, moderate-severe turbulence. Load
%   factor commands are capped by the lift available (0.85 CLmax qS/W), so
%   speed bleed shows up as reduced g rather than a stall.

d2r = pi/180;  g0 = 9.80665;
P.V0 = 250;  P.h0 = 5000;  P.hdg0 = 90;  P.model = '6dof';
P.aircraft = 'F16';
P.jerkGps = 15;  P.pddMax = 3000*d2r;  P.wind = [];  P.dCl = [];  P.vibe = [];
V0 = 250;
ph = struct('name',{},'cmd',{},'done',{},'maneuver',{});
acM = [9295 27.87 1.89];                                   % F-16 [mass kg, wing m^2, CLmax]
nAvail = @(s) 0.85*acM(3)*0.5*rhoStd(s.h)*s.V^2*acM(2)/(acM(1)*g0);   % lift-limited load factor
% aileron sign giving a RIGHT roll (F-16 tables: positive aileron rolls left)
sA = -1;

switch name
    case 'baseline_cruise'
        P.description = 'Straight and level, 250 m/s at 5000 m (reference)';
        ph(end+1) = cruise('cruise', 60, V0, true);

    case 'turn'
        P.description = 'Max-performance level turn: 3 g/s onset to 9 g (lift-limited as speed bleeds), 22 s';
        Ton = 8*pi/(2*3);  hold = 22;
        ph(end+1) = cruise('cruise', 12, V0, false);
        ph(end+1) = mk('turn9g', @(tau,s,s0) turnCmd(tau,s,s0,Ton,hold,9,+1,V0,nAvail), @(tau,s,s0) tau >= 2*Ton + hold, true);
        ph(end+1) = cruise('cruise', 60, V0, false);

    case 'rollin'
        P.description = 'Turn entries/exits: max-rate roll (270 deg/s requested) to 80 deg bank, ~6-7 g level turns, L/R x4';
        ph(end+1) = cruise('cruise', 12, V0, false);
        for k = 1:4
            dir = (-1)^(k+1);
            ph(end+1) = rollTo(sprintf('rollin%d',k), 80*dir, 270, 4.0, V0, @(s,s0,tau) min(blendN(nLevel(s,s0.h),s0,tau,0.5), 6.1), true); %#ok<AGROW>
            ph(end+1) = unloadPh(sprintf('unload%d',k), 1.5, 0.8, V0); %#ok<AGROW>
            ph(end+1) = rollTo(sprintf('rollout%d',k), 0, 270, 4.5, V0, @(s,s0,tau) min(nLevel(s,s0.h), 2.5), true); %#ok<AGROW>
        end
        ph(end+1) = cruise('cruise', 60, V0, false);

    case 'loop_pushover'
        P.description = 'Loop at up to 7 g (lift-limited), then 4 g pull-up to 30 deg and -2 g push-over';
        ph(end+1) = cruise('cruise', 12, V0, false);
        ph(end+1) = mk('loop', @(tau,s,s0) o6(0, min(1 + 6*rc(tau,1.57), nAvail(s)), vHold(s,V0), false, 0, []), ...
                       @(tau,s,s0) s.qint >= 2*pi - 0.35, true);
        ph(end+1) = cruise('recover', 4, V0, true);
        ph(end+1) = mk('pullup', @(tau,s,s0) o6(bankP(s,0,0), blendN(4,s0,tau,0.8), vHold(s,V0), false, 0, []), ...
                       @(tau,s,s0) s.gamma >= 30*d2r, true);
        ph(end+1) = mk('pushover', @(tau,s,s0) o6(bankP(s,0,0), blendN(-2,s0,tau,0.8), vHold(s,V0), false, 0, []), ...
                       @(tau,s,s0) s.gamma <= -20*d2r, true);
        ph(end+1) = cruise('recover', 5, V0, true);
        ph(end+1) = cruise('cruise', 60, V0, false);

    case 'break_gonset'
        P.description = 'Break turns: 200 deg/s roll to 84 deg, 12 g/s onset to 9 g, 5 s hold, R then L';
        ph(end+1) = cruise('cruise', 12, V0, false);
        for dir = [1 -1]
            phiB = acos(1/9)/d2r*dir;
            ph(end+1) = rollTo('roll', phiB, 200, 0, V0, @(s,s0,tau) blendN(1,s0,tau,0.5), true); %#ok<AGROW>
            ph(end+1) = mk('pull9g', @(tau,s,s0) o6(bankP(s,phiB*d2r,0), min(1 + 8*rc(tau, 8*pi/(2*12)), nAvail(s)), ...
                           vHold(s,V0), false, 0, []), @(tau,s,s0) tau >= 6.05, true); %#ok<AGROW>
            ph(end+1) = unloadPh('unload', 1.5, 0.8, V0); %#ok<AGROW>
            ph(end+1) = rollTo('rollout', 0, 200, 6, V0, @(s,s0,tau) min(nLevel(s,s0.h), 2.5), true); %#ok<AGROW>
        end
        ph(end+1) = cruise('cruise', 60, V0, false);

    case 'jinking'
        P.description = 'Evasive jinks: 6 g level-turn reversals every ~2.7 s (300 deg/s), then +5/-2 g pitch jinks';
        ph(end+1) = cruise('cruise', 12, V0, false);
        phiJ = acos(1/6)/d2r;
        ph(end+1) = rollTo('jink1', phiJ, 300, 1.8, V0, [], true);
        for k = 2:6
            ph(end+1) = rollTo(sprintf('jink%d',k), phiJ*(-1)^(k+1), 300, 1.8, V0, [], true); %#ok<AGROW>
        end
        ph(end+1) = unloadPh('unload', 1.5, 0.6, V0);
        ph(end+1) = rollTo('level', 0, 300, 2, V0, @(s,s0,tau) min(nLevel(s,s0.h), 2.5), true);
        for k = 1:3
            ph(end+1) = mk('pull5', @(tau,s,s0) o6(bankP(s,0,0), min(blendN(5,s0,tau,0.4), nAvail(s)), vHold(s,V0), false, 0, []), ...
                           @(tau,s,s0) tau >= 1.2, true); %#ok<AGROW>
            ph(end+1) = mk('push2', @(tau,s,s0) o6(bankP(s,0,0), blendN(-2,s0,tau,0.5), vHold(s,V0), false, 0, []), ...
                           @(tau,s,s0) tau >= 1.2, true); %#ok<AGROW>
        end
        ph(end+1) = cruise('recover', 4, V0, true);
        ph(end+1) = cruise('cruise', 60, V0, false);

    case 'rolls'
        P.description = 'Three consecutive full-aileron rolls, then a 2-turn max-rate rudder-assisted (snap-style) roll';
        ph(end+1) = cruise('cruise', 15, V0, false);
        ph(end+1) = rollBy('aileron3', 1080, 300, 0.25, V0, true);
        ph(end+1) = cruise('recover', 8, V0, true);
        % full aileron + half rudder, elevator left to the autopilot (n = 1): an
        % open-loop snap roll at 230 m/s over-stressed the jet (12.5 g)
        ph(end+1) = mk('snap', @(tau,s,s0) o6(0, 1, vHold(s,V0), false, 0, [NaN, sA*21.5*d2r, -15*d2r, NaN]), ...
                       @(tau,s,s0) s.pint >= 4*pi - 0.6 || tau >= 5, true);
        ph(end+1) = cruise('recover', 6, V0, true);
        ph(end+1) = cruise('cruise', 60, V0, false);

    case 'splitS_immelmann'
        P.description = 'Split-S (half roll, 6 g pull through) then Immelmann (5 g half loop, half roll)';
        ph(end+1) = cruise('cruise', 12, V0, false);
        ph(end+1) = rollBy('halfroll', 180, 180, 0.2, V0, true);
        ph(end+1) = mk('splitS', @(tau,s,s0) o6(0, min(blendN(6,s0,tau,0.8), nAvail(s)), vHold(s,V0), false, 0, []), ...
                       @(tau,s,s0) s.qint >= pi - 0.1, true);
        ph(end+1) = cruise('recover', 6, V0, true);
        ph(end+1) = mk('immelmann', @(tau,s,s0) o6(0, min(blendN(5,s0,tau,0.8), nAvail(s)), vHold(s,V0), false, 0, []), ...
                       @(tau,s,s0) s.qint >= pi - 0.1, true);
        ph(end+1) = rollBy('halfroll', 180, 180, 0.2, V0, true);
        ph(end+1) = cruise('recover', 6, V0, true);
        ph(end+1) = cruise('cruise', 60, V0, false);

    case 'departure_spin'
        P.description = 'Idle 5 g turn to bleed to 140 m/s, full aft stick + pro-spin rudder/aileron 5 s, anti-spin recovery, pull-out';
        ph(end+1) = cruise('cruise', 12, V0, false);
        % speed bled in an idle level turn (an idle pull-up zoomed to 7 km and
        % tail-slid to 16 m/s, far outside the F-16 data's 91 m/s lower bound)
        ph(end+1) = rollTo('bank78', 78.5, 90, 0, V0, [], true);
        ph(end+1) = mk('bleed', @(tau,s,s0) o6(bankP(s,78.5*d2r,0), min(blendN(nLevel(s,s0.h),s0,tau,0.8), nAvail(s)), -20, false, 0, [NaN NaN NaN 0]), ...
                       @(tau,s,s0) s.V <= 140, true);
        ph(end+1) = mk('departure', @(tau,s,s0) o6(0, 1, 0, false, 0, [-25*d2r, sA*20*d2r, -30*d2r, 0]), ...
                       @(tau,s,s0) tau >= 5, true);
        ph(end+1) = mk('antispin', @(tau,s,s0) o6(0, 1, 0, false, 0, [8*d2r, 0, 30*d2r, 0]), ...
                       @(tau,s,s0) (tau >= 1 && abs(s.w(3)) < 15*d2r && abs(s.w(1)) < 30*d2r) || tau >= 4, true);
        ph(end+1) = mk('wingslevel', @(tau,s,s0) o6(bankP(s,0,0), 1, 0, false, 0, [NaN NaN NaN 1]), ...
                       @(tau,s,s0) abs(s.phi) < 10*d2r || tau >= 3, true);
        ph(end+1) = mk('pullout', @(tau,s,s0) o6(bankP(s,0,0), min(blendN(5,s0,tau,0.8), nAvail(s)), vHold(s,V0), false, 0, []), ...
                       @(tau,s,s0) tau >= 1 && s.gamma >= 0, true);
        ph(end+1) = cruise('recover', 4, V0, true);
        ph(end+1) = cruise('cruise', 60, V0, false);

    case 'catapult_trap'
        P.description = 'Deck, 4.5 g catapult to 85 m/s, climb/4 g turn/approach, arrested landing ~4 g';
        P.model = 'pointmass';  P.V0 = 0;  P.h0 = 20;  P.jerkGps = 60;
        ph(end+1) = mk('deck', @(tau,s,s0) o6(0,1,0,true,0,[]), @(tau,s,s0) tau >= 10, false);
        ph(end+1) = mk('catapult', @(tau,s,s0) o6(0,1, 45*rc(tau,0.15)*(1 - rc(tau-1.75,0.15)), true, 0, []), ...
                       @(tau,s,s0) tau >= 2.05, true);
        ph(end+1) = mk('climb', @(tau,s,s0) o6(bankP(s,0,0), nGamma(s, 12*d2r*rc(tau,1.5), 2), vHold(s,100), false, 0, []), ...
                       @(tau,s,s0) s.h >= 150, true);
        ph(end+1) = mk('level', @(tau,s,s0) o6(bankP(s,0,0), nLevel(s,150), vHold(s,100), false, 0, []), ...
                       @(tau,s,s0) tau >= 3, true);
        ph(end+1) = rollTo('turn4g', acos(1/4)/d2r, 90, 0, 100, @(s,s0,tau) nLevel(s,150), true);
        ph(end).done = @(tau,s,s0) wrapPi(s.psi - s0.psi) >= 165*d2r;
        ph(end+1) = rollTo('rollout', 0, 90, 2, 100, @(s,s0,tau) nLevel(s,150), true);
        ph(end+1) = mk('approach', @(tau,s,s0) o6(bankP(s,0,0), nGamma(s,-8*d2r*rc(tau,3),1), vHold(s,75), false, 0, []), ...
                       @(tau,s,s0) s.h <= 22, true);
        ph(end+1) = mk('touchdown', @(tau,s,s0) o6(0, nGamma(s,0,8), 0, false, 0, []), ...
                       @(tau,s,s0) s.gamma >= -0.1*d2r || tau >= 0.6, true);
        ph(end+1) = mk('trap', @(tau,s,s0) o6(0,1, -pi*s0.V/(2*2.75)*sin(pi*min(tau,2.75)/2.75), true, 0, []), ...
                       @(tau,s,s0) tau >= 2.75, true);
        ph(end+1) = mk('deck', @(tau,s,s0) o6(0,1,0,true,0,[]), @(tau,s,s0) false, false);

    case 'hard_landing'
        P.description = 'Approach at 75 m/s on a 3 deg path, hard touchdown (~3 g spike), braking rollout';
        P.model = 'pointmass';  P.V0 = 75;  P.h0 = 120;  P.jerkGps = 60;
        ph(end+1) = cruise('approachlevel', 10, 75, false);
        ph(end+1) = mk('glide', @(tau,s,s0) o6(bankP(s,0,0), nGamma(s,-3*d2r*rc(tau,2),2), vHold(s,75), false, 0, []), ...
                       @(tau,s,s0) s.h <= 20.5, false);
        ph(end+1) = mk('touchdown', @(tau,s,s0) o6(0, nGamma(s,0,8), 0, false, 0, []), ...
                       @(tau,s,s0) s.gamma >= -0.05*d2r || tau >= 0.8, true);
        ph(end+1) = mk('rollout', @(tau,s,s0) o6(0,1,-3.5,true,0,[]), @(tau,s,s0) s.V <= 15, true);
        ph(end+1) = mk('taxi', @(tau,s,s0) o6(0,1,0,true,0,[]), @(tau,s,s0) false, false);

    case 'turbulence'
        P.description = 'Moderate-severe turbulence (gusts 3/4/6 m/s N/E/D) 12-48 s, wake-vortex roll upset at 30 s, vibration';
        g1 = smoothNoise(11, 0.1, 2, 14);  g2 = smoothNoise(12, 0.1, 2, 14);  g3 = smoothNoise(13, 0.1, 2.5, 14);
        win = @(t) rc(t-12, 1).*(1 - rc(t-47, 1));
        P.wind = @(t) win(t)*[3*g1(t); 4*g2(t); 6*g3(t)];
        P.dCl  = @(t) 0.04*sin(pi*min(max(t-30,0),0.8)/0.8)^2;
        P.vibe = @(t) win(t(:)).*g0.*[0.08*sin(2*pi*31*t(:)), 0.05*sin(2*pi*43*t(:)+0.5), ...
                                    0.25*sin(2*pi*22*t(:)) + 0.15*sin(2*pi*67*t(:)+1)];
        ph(end+1) = cruise('cruise', 12, V0, false);
        ph(end+1) = cruise('turbulence', 36, V0, true);
        ph(end+1) = cruise('cruise', 60, V0, false);

    otherwise
        if ~strncmp(name, 'custom', 6)
            error('fighterPhases: unknown profile "%s" (see fighterProfileList).', name);
        end
        % user-defined manoeuvre: opt.customPhases is an N x 1 cell array of
        % building blocks (see customManeuverBlocks.m); a 12 s straight-and-
        % level lead-in (shared with the baseline, for the receiver's warm
        % start) and a trailing cruise to the end of the record are added
        if nargin < 2 || ~isfield(opt,'customPhases') || isempty(opt.customPhases)
            error('fighterPhases: custom profile "%s" needs opt.customPhases.', name);
        end
        P.description = 'Custom manoeuvre';
        if isfield(opt,'customDescription') && ~isempty(opt.customDescription)
            P.description = opt.customDescription;
        end
        ph(end+1) = cruise('cruise', 12, V0, false);
        for k = 1:numel(opt.customPhases)
            b = opt.customPhases{k};  typ = lower(b{1});  a = [b{2:end}];
            switch typ
                case 'cruise'                       % {'cruise', T}
                    ph(end+1) = cruise('cruise', a(1), V0, true); %#ok<AGROW>
                case 'turn'                         % {'turn', n, hold, onset, dir}
                    nPk = a(1);  hold = a(2);  onset = a(3);  dir = sign(a(4));  if dir == 0, dir = 1; end
                    Ton = (nPk - 1)*pi/(2*onset);
                    % load factor capped at the requested peak (the level-hold of
                    % turnCmd overshoots during a fast bank onset), then roll out
                    % to wings level so that the next block starts level
                    ph(end+1) = mk(sprintf('turn%gg', nPk), @(tau,s,s0) turnCmdCapped(tau,s,s0,Ton,hold,nPk,dir,V0,nAvail), ...
                                   @(tau,s,s0) tau >= 2*Ton + hold, true); %#ok<AGROW>
                    ph(end+1) = rollTo('rollout', 0, 200, 1.0, V0, @(s,s0,tau) min(nLevel(s,s0.h), 2.5), true); %#ok<AGROW>
                case 'bank'                         % {'bank', phi, rate, hold}
                    % level-hold load factor, capped near what level flight at
                    % the target bank needs (the uncapped level-hold demanded
                    % ~9 g when rolling level after a climb)
                    nCapB = max(2.5, 1.15/cos(min(abs(a(1)), 85)*d2r));
                    ph(end+1) = rollTo(sprintf('bank%g', a(1)), a(1), a(2), a(3), V0, ...
                                       @(s,s0,tau) min(blendN(nLevel(s,s0.h),s0,tau,0.5), nCapB), true); %#ok<AGROW>
                case 'roll'                         % {'roll', angle, rate}
                    ph(end+1) = rollBy(sprintf('roll%g', a(1)), a(1), a(2), 0.25, V0, true); %#ok<AGROW>
                case 'pull'                         % {'pull', n, T}  (n < 0: push)
                    nT = a(1);  Tp = a(2);
                    ph(end+1) = mk(sprintf('pull%gg', nT), @(tau,s,s0) o6(bankP(s,0,0), min(blendN(nT,s0,tau,0.8), nAvail(s)), ...
                                   vHold(s,V0), false, 0, []), @(tau,s,s0) tau >= Tp, true); %#ok<AGROW>
                case 'loop'                         % {'loop', n}
                    nL = a(1);
                    ph(end+1) = mk('loop', @(tau,s,s0) o6(0, min(1 + (nL - 1)*rc(tau,1.57), nAvail(s)), vHold(s,V0), false, 0, []), ...
                                   @(tau,s,s0) s.qint >= 2*pi - 0.35, true); %#ok<AGROW>
                    ph(end+1) = cruise('recover', 4, V0, true); %#ok<AGROW>
                case 'unload'                       % {'unload', n, T}
                    ph(end+1) = unloadPh('unload', a(1), a(2), V0); %#ok<AGROW>
                otherwise
                    error('fighterPhases: unknown custom block "%s" (see customManeuverBlocks).', typ);
            end
        end
        ph(end+1) = cruise('recover', 4, V0, true);
        ph(end+1) = cruise('cruise', 60, V0, false);
end
P.phases = ph;
end

%% ---------------------------------------------------------------- phases
function ph = mk(name, cmd, done, man)
ph = struct('name', name, 'cmd', cmd, 'done', done, 'maneuver', man);
end

function ph = cruise(name, dur, Vref, man)
% gentle: load factor capped at 2.5 g (recovering from a steep bank under
% level-hold otherwise demanded ~5 g while rolling back)
ph = mk(name, @(tau,s,s0) o6(bankP(s,0,0), min(blendN(nLevel(s,s0.h),s0,tau,0.8), 2.5), vHold(s,Vref), false, 0, []), ...
        @(tau,s,s0) tau >= dur, man);
end

function ph = rollTo(name, phiTdeg, pMaxDeg, holdDur, Vref, nFn, man)
% roll to bank phiT with a raised-cosine bank history (peak rate pMax), then hold
d2r = pi/180;  phiT = phiTdeg*d2r;  pMax = pMaxDeg*d2r;
if isempty(nFn), nFn = @(s,s0,tau) blendN(nLevel(s,s0.h),s0,tau,0.5); end
% duration: peak rate pMax, peak roll acceleration (pi/Tr)^2*dphi/2 within
% aMax = 10 rad/s^2 and peak roll jerk (pi/Tr)^3*dphi/2 within jMax = 15
% rad/s^3 - the binding limit of the F-16's 80 deg/s aileron actuators
% (a full aileron reversal takes 0.54 s)
aMax = 10;  jMax = 15;
Tr = @(s0) max([abs(phiT - s0.phi)*pi/(2*pMax), pi*sqrt(abs(phiT - s0.phi)/(2*aMax)), ...
                pi*(abs(phiT - s0.phi)/(2*jMax))^(1/3), 0.05]);
ph = mk(name, @(tau,s,s0) o6(rollRef(tau,s,s0,phiT,Tr(s0)), nFn(s,s0,tau), vHold(s,Vref), false, 0, []), ...
        @(tau,s,s0) tau >= Tr(s0) + holdDur, man);
end

function ph = unloadPh(name, nT, Tb, Vref)
% hold the current bank and unload to nT over Tb before rolling out (rolling
% out of a steep bank while still loaded drives alpha past the limit)
ph = mk(name, @(tau,s,s0) o6(bankP(s,s0.phi,0), blendN(nT,s0,tau,Tb), vHold(s,Vref), false, 0, []), ...
        @(tau,s,s0) tau >= Tb, true);
end

function p = rollRef(tau, s, s0, phiT, Tr)
dphi = phiT - s0.phi;
phiC = s0.phi + dphi*rc(tau, Tr);
phiD = 0;
if tau > 0 && tau < Tr, phiD = dphi*pi/(2*Tr)*sin(pi*tau/Tr); end
p = bankP(s, phiC, phiD);
end

function ph = rollBy(name, dphiDeg, pMaxDeg, Ta, Vref, man)
% roll through dphi (accumulated body roll angle) at up to pMax with n = 1:
% raised-cosine rate onset over Ta, then a rate limited by the remaining
% angle (sqrt(2 a rem), a = 8 rad/s^2) so the roll ends where it should even
% when the aircraft cannot reach pMax (a time-scheduled roll ended 93 deg short)
d2r = pi/180;  sg = sign(dphiDeg);  pMax = pMaxDeg*d2r;  dphi = abs(dphiDeg)*d2r;
rem = @(s) dphi - sg*s.pint;
ph = mk(name, @(tau,s,s0) o6(sg*min(pMax*rc(tau,Ta), sqrt(2*8*max(rem(s),0))), 1, vHold(s,Vref), false, 0, []), ...
        @(tau,s,s0) rem(s) <= 0.02 || tau >= dphi/pMax*3 + 2, man);
end

function varargout = turnCmd(tau, s, s0, Ton, hold, nPk, dir, Vref, nAvail)
% level turn: load factor schedule capped by the available lift; the bank
% follows the achievable load factor (cos(phi) = 1/n) so the turn stays
% level and bleeds speed towards corner speed like a real max-rate turn
nCap = max(nAvail(s), 1.05);
nT  = @(t) min(1 + (nPk - 1)*(rc(t,Ton) - rc(t - Ton - hold, Ton)), nCap);
phc = @(t) dir*acos(1/nT(t));
h = 1e-3;
p = bankP(s, phc(tau), (phc(tau+h) - phc(tau-h))/(2*h));
out = {p, min(nLevel(s, s0.h), nCap), vHold(s, Vref), false, 0, []};
varargout = out(1:nargout);
end

function varargout = turnCmdCapped(tau, s, s0, Ton, hold, nPk, dir, Vref, nAvail)
% turnCmd with the load factor limited to 1.05 x the requested peak (custom turns)
[p, n, Vd, gnd, nDir, ovr] = turnCmd(tau, s, s0, Ton, hold, nPk, dir, Vref, nAvail);
out = {p, min(n, 1.05*nPk), Vd, gnd, nDir, ovr};
varargout = out(1:nargout);
end

%% --------------------------------------------------------------- helpers
function varargout = o6(p, n, Vd, gnd, nDir, ovr)
out = {p, n, Vd, gnd, nDir, ovr};
varargout = out(1:nargout);
end

function y = rc(tau, T)
% raised-cosine step 0 -> 1 over [0, T]
y = 0.5*(1 - cos(pi*min(max(tau./T, 0), 1)));
end

function a = wrapPi(a)
a = mod(a + pi, 2*pi) - pi;
end

function p = bankP(s, phiC, phiD)
p = phiD + 3*wrapPi(phiC - s.phi);
end

function n = blendN(nTarget, s0, tau, Tb)
% move from the load factor at the phase start to nTarget over Tb (bounded onset)
n = s0.nAct + (nTarget - s0.nAct)*rc(tau, Tb);
end

function n = nLevel(s, hRef)
g0 = 9.80665;  V = max(s.V, 30);
gc = asin(max(min(0.25*(hRef - s.h)/V, 0.25), -0.25));
n  = nGammaCore(s, gc, 2.0, V, g0);
end

function n = nGamma(s, gc, K)
n = nGammaCore(s, gc, K, max(s.V, 30), 9.80665);
end

function n = nGammaCore(s, gc, K, V, g0)
% never flip the sign of the demand: past 90 deg of bank (e.g. a reversal
% overshoot) "level" would mean a negative-g push, which no pilot flies
% mid-reversal; intentionally inverted phases set n explicitly
% flight-path correction (in g) limited to +-0.7 g when banked beyond 30 deg:
% in a steep bank the vertical response to n is weak (n cos(phi)), and an
% unlimited correction for a small climb in a 6 g reversal demanded -3 g
c = max(cos(s.phi), 0.2);
corr = V*K*(gc - s.gamma)/g0;
if abs(s.phi) > 30*pi/180, corr = max(min(corr, 0.7), -0.7); end
n = (cos(s.gamma) + corr)/c;
n = max(min(n, 9.5), -3);
end

function Vd = vHold(s, Vref)
Vd = max(min(0.4*(Vref - s.V), 8), -8);
end

function r = rhoStd(h)
T = 288.15 - 0.0065*min(max(h,0),11000);
r = 1.225*(T/288.15)^4.2559;
end

function g = smoothNoise(seed, fmin, fmax, nComp)
% unit-variance smooth random signal: sum of sines, log-spaced frequencies
st = rng;  rng(seed);
f  = logspace(log10(fmin), log10(fmax), nComp);
ph = 2*pi*rand(1, nComp);
rng(st);
a  = sqrt(2/nComp);
g  = @(t) sum(a*sin(2*pi*f(:)*t(:)' + ph(:)), 1)';
end
