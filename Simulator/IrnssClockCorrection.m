function satClkCorr = IrnssClockCorrection(transmitTime, eph)
% Satellite clock correction for the NavIC/IRNSS L5 SPS signal.
%
% Mirrors the formula used by NavIC_Software_Receiver/include/satpos.m so
% the simulated satellite clock error matches what the receiver will
% recompute from the same broadcast af0/af1/af2/TGD: the broadcast TGD is
% referenced to the S-band signal, so it is scaled by
% JS_gamma = (F_S/F_L5)^2 for use on L5, per the IRNSS ICD.
%
%--- Find time difference
dt = check_t(transmitTime - eph.toc);

F_S  = 2492.028e6;
F_L5 = 1176.45e6;
JS_gamma = (F_S / F_L5)^2;

satClkCorr = eph.af0 + eph.af1 * dt + eph.af2 * dt * dt - JS_gamma * eph.tgd;
