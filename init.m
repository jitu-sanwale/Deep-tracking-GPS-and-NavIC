%INIT  GPS L1 C/A and NavIC L5 SPS: IF signal simulator and software receiver.
%
%   Menu:
%     1) generate an IF dataset         (Simulator\init.m)
%     2) process an IF dataset          (Receiver\init.m)
%     3) both: generate, then process
%
%   Simulator: truth trajectories for predefined or custom aircraft
%   manoeuvres (nonlinear F-16 6-DOF model) or a user CSV trajectory, a
%   navigation-grade IMU record, and GPS L1 C/A and NavIC L5 SPS IF files
%   (5 MHz complex, 8-bit I/Q; 600 MB per 60 s file and band).
%   Receiver: classical scalar tracking (CT), vector tracking (VT) and
%   INS-aided ultratight tracking (UT) in GPS-only, NavIC-only and combined
%   GPS + NavIC configurations, with comparison against the truth trajectory.
%
%   Both parts also run on their own from Simulator\init.m and
%   Receiver\init.m; see their headers for non-interactive use (SIM and RX
%   structs). Requirements: MATLAB R2023a or later with the Signal
%   Processing, Communications and Curve Fitting Toolboxes.
%
%   The receiver builds on GPSSDR_vt (B. Xu and L.-T. Hsu, GPS Solutions 23,
%   46, 2019) and the simulator on GPS_L1CA_signal_generator, the open-source
%   IF signal simulator of Y. Li (Beijing Information Science and Technology
%   University) as modified by J. Sun; both include routines from SoftGNSS
%   (K. Borre, D. Akos, D. Plausinaitis et al.).
%   The F-16 aerodynamic data are from F16Model.jl (MIT licence, see
%   Simulator\F16Aircraft). The broadcast ephemerides
%   (Simulator\IITK00IND_R_20250011000_01H_MN.rnx) are data of the
%   International GNSS Service (IGS), https://igs.org.
%   The original notices are kept in the file headers. Licence: GPL-3.0-or-later.

pkg = fileparts(mfilename('fullpath'));
fprintf(['\nGPS / NavIC IF simulator and software receiver (%s)\n' ...
         '  1) generate an IF dataset (simulator)\n' ...
         '  2) process an IF dataset (receiver)\n' ...
         '  3) both: generate, then process\n'], pkg);
choice = strtrim(input('Select [3]: ', 's'));
if isempty(choice), choice = '3'; end
switch choice
    case '1'
        clear SIM;  run(fullfile(pkg, 'Simulator', 'init.m'));
    case '2'
        clear RX;   run(fullfile(pkg, 'Receiver', 'init.m'));
    case '3'
        clear SIM RX;
        run(fullfile(pkg, 'Simulator', 'init.m'));
        % the receiver reads the data folder the simulator has just written
        RX = struct('dataDir', SIMDATADIR);
        run(fullfile(pkg, 'Receiver', 'init.m'));
    otherwise
        error('init: choose 1, 2 or 3.');
end
