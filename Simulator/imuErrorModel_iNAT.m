function imu = imuErrorModel_iNAT(tF, fB, wB, fsIMU, nSub, accelRangeG, seed)
%IMUERRORMODEL_INAT iMAR iNAT-FSLG-01 sensor error model applied to truth.
%   imu = imuErrorModel_iNAT(tF, fB, wB, fsIMU, nSub, accelRangeG, seed)
%
%   tF      - internal truth grid [s], nOut*nSub+1 points, spacing 1/(fsIMU*nSub)
%   fB, wB  - true specific force [m/s^2] and angular rate [rad/s] (body)
%             on tF, N x 3
%   fsIMU   - output rate [Hz];  nSub - even number of Simpson sub-intervals
%   accelRangeG - accelerometer range [g];  seed - rng seed
%
%   Every term of the iNAT-FSLG-01/T datasheet (rev 2.29, rms values, each
%   treated as a 1-sigma per-axis draw):
%                             Gyroscopes            Accelerometers
%       range                 +/-600 deg/s          +/-accelRangeG (DS: 20 g)
%       white noise           ARW 0.0045 deg/rt-h   40 ug/rt-Hz
%       turn-on bias          0.05 deg/h            300 ug
%       bias instability      0.01 deg/h            15 ug      (1st-order
%                             Gauss-Markov, tau = 300 s)
%       scale factor          150 ppm (0.015 %)     300 ppm (0.03 %)
%       non-linearity         100 ppm of FS (0.01 %) 50 ug/g^2
%       axis misalignment     0.15 mrad             0.25 mrad
%       resolution            0.05 urad/LSB         0.1 ug/LSB
%   Order: true -> misalignment & scale factor -> non-linearity -> range
%   clipping (on the raw grid, like the sensor itself) -> integrate (Simpson)
%   -> bias + instability + white noise -> quantisation with residual carry.
%   Outputs f, w are interval averages over (t-dt, t].

g0 = 9.80665;  d2r = pi/180;
dtIMU = 1/fsIMU;
nOut  = (numel(tF) - 1)/nSub;

%% Datasheet error parameters
ds.gyroRange    = 600*d2r;              % rad/s
ds.gyroARW      = 0.0045*d2r/60;        % rad/sqrt(s)
ds.gyroBias     = 0.05*d2r/3600;        % rad/s
ds.gyroBI       = 0.01*d2r/3600;        % rad/s
ds.gyroSF       = 150e-6;
ds.gyroNL       = 100e-6;               % of full scale
ds.gyroMis      = 0.15e-3;              % rad
ds.gyroLSB      = 0.05e-6;              % rad per LSB of dtheta
ds.accRange     = accelRangeG*g0;       % m/s^2
ds.accVRW       = 40e-6*g0;             % m/s^2/sqrt(Hz) = m/s^1.5
ds.accBias      = 300e-6*g0;            % m/s^2
ds.accBI        = 15e-6*g0;             % m/s^2
ds.accSF        = 300e-6;
ds.accNL        = 50e-6/g0;             % (m/s^2) per (m/s^2)^2
ds.accMis       = 0.25e-3;              % rad
ds.accLSB       = 0.1e-6*g0;            % m/s^2 per LSB
ds.tauBI        = 300;                  % s, Gauss-Markov correlation time

%% Error realisation
rng(seed);
offd = @(s) s*(randn(3).*~eye(3));
err.Mg   = diag(ds.gyroSF*randn(3,1)) + offd(ds.gyroMis);
err.Ma   = diag(ds.accSF*randn(3,1))  + offd(ds.accMis);
err.kNLg = ds.gyroNL*randn(3,1);
err.kNLa = ds.accNL*randn(3,1);
err.bg0  = ds.gyroBias*randn(3,1);
err.ba0  = ds.accBias*randn(3,1);

%% Deterministic sensor errors and clipping on the raw grid
fRaw = fB*(eye(3) + err.Ma)' + (fB.^2).*err.kNLa';
wRaw = wB*(eye(3) + err.Mg)' + (wB.*abs(wB)/ds.gyroRange).*err.kNLg';
clipA = abs(fRaw) > ds.accRange;   clipG = abs(wRaw) > ds.gyroRange;
fRaw = max(min(fRaw, ds.accRange), -ds.accRange);
wRaw = max(min(wRaw, ds.gyroRange), -ds.gyroRange);

%% Integrate to increments (Simpson) and average
wS = [1, repmat([4 2], 1, nSub/2-1), 4, 1]' * (dtIMU/nSub)/3;
dvT = zeros(nOut,3); dthT = zeros(nOut,3); dv = zeros(nOut,3); dth = zeros(nOut,3);
clipOut = false(nOut,1);
for j = 1:nOut
    idx = (j-1)*nSub + (1:nSub+1);
    dvT(j,:)  = wS'*fB(idx,:);   dthT(j,:) = wS'*wB(idx,:);
    dv(j,:)   = wS'*fRaw(idx,:); dth(j,:)  = wS'*wRaw(idx,:);
    clipOut(j) = any(any(clipA(idx,:) | clipG(idx,:)));
end
tOut = tF(nSub+1:nSub:end);

%% Stochastic errors: bias, Gauss-Markov instability, white noise
phiGM = exp(-dtIMU/ds.tauBI);
qGM   = sqrt(1 - phiGM^2);
bgGM = zeros(nOut,3); baGM = zeros(nOut,3);
bgGM(1,:) = ds.gyroBI*randn(1,3);  baGM(1,:) = ds.accBI*randn(1,3);
for j = 2:nOut
    bgGM(j,:) = phiGM*bgGM(j-1,:) + qGM*ds.gyroBI*randn(1,3);
    baGM(j,:) = phiGM*baGM(j-1,:) + qGM*ds.accBI*randn(1,3);
end
dth = dth + (err.bg0' + bgGM)*dtIMU + ds.gyroARW*sqrt(dtIMU)*randn(nOut,3);
dv  = dv  + (err.ba0' + baGM)*dtIMU + ds.accVRW*sqrt(dtIMU)*randn(nOut,3);

%% Quantisation with residual carry
qth = ds.gyroLSB;  qv = ds.accLSB*dtIMU;
resTh = zeros(1,3); resV = zeros(1,3);
for j = 1:nOut
    x = dth(j,:) + resTh;  dth(j,:) = round(x/qth)*qth;  resTh = x - dth(j,:);
    x = dv(j,:)  + resV;   dv(j,:)  = round(x/qv)*qv;    resV  = x - dv(j,:);
end

imu.t         = tOut;
imu.fs        = fsIMU;
imu.f         = dv/dtIMU;        imu.w      = dth/dtIMU;
imu.dv        = dv;              imu.dtheta = dth;
imu.fTrue     = dvT/dtIMU;       imu.wTrue  = dthT/dtIMU;
imu.clipped   = clipOut;
imu.biasG     = err.bg0' + bgGM; imu.biasA  = err.ba0' + baGM;
imu.err       = err;             imu.datasheet = ds;
imu.accelRangeG = accelRangeG;   imu.seed = seed;
imu.frame     = 'body x-fwd y-right z-down; f, w are interval averages over (t-dt,t]';
end
