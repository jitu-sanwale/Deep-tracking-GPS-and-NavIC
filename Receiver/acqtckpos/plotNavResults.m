function plotNavResults(navCT, navVT, CN0_CT, CN0_VT, svList, file, outDir)
%Purpose:
%   Plot conventional (CT) and vector tracking (VT) navigation results over
%   the whole data file and save figures (PNG) and a summary CSV in outDir.
if ~exist(outDir,'dir'), mkdir(outDir); end
hasVT = ~isempty(navVT);

tCT = (1:size(navCT.usrPosENU,1))'*1e-3;                      % s from start of positioning
if hasVT
    tVT = navVT.tData(:) - file.skip*1e-3;                 % VT solutions every track.utUpdateMs
end
cCT = [0.85 0.33 0.10]; cVT = [0 0.45 0.74];
lbl = {'East','North','Up'};

%% 1) ENU position vs time
f1 = figure('Name','Position ENU','Color','w','Position',[50 50 1000 750]);
for k = 1:3
    subplot(3,1,k); hold on; grid on;
    plot(tCT, navCT.usrPosENU(:,k), 'Color', cCT);
    if hasVT, plot(tVT, navVT.usrPosENU(:,k), 'Color', cVT); end
    ylabel([lbl{k} ' (m)']);
    if k == 1, title('Position relative to reference (ENU)'); legend(legendNames(hasVT),'Location','best'); end
end
xlabel('Time since start of positioning (s)');
saveFig(f1, fullfile(outDir,'01_position_ENU.png'));

%% 2) Horizontal trajectory
f2 = figure('Name','Trajectory','Color','w','Position',[80 80 800 750]); hold on; grid on; axis equal;
plot(navCT.usrPosENU(:,1), navCT.usrPosENU(:,2), 'Color', cCT);
if hasVT, plot(navVT.usrPosENU(:,1), navVT.usrPosENU(:,2), 'Color', cVT); end
plot(navCT.usrPosENU(1,1), navCT.usrPosENU(1,2), 'k^', 'MarkerFaceColor','k');
xlabel('East (m)'); ylabel('North (m)'); title('Horizontal trajectory');
legend([legendNames(hasVT) {'Start'}],'Location','best');
saveFig(f2, fullfile(outDir,'02_trajectory_EN.png'));

%% 3) Velocity
f3 = figure('Name','Velocity','Color','w','Position',[110 110 1000 750]);
for k = 1:3
    subplot(4,1,k); hold on; grid on;
    plot(tCT, navCT.usrVelENU(:,k), 'Color', cCT);
    if hasVT, plot(tVT, navVT.usrVelENU(:,k), 'Color', cVT); end
    ylabel(['V_' lbl{k}(1) ' (m/s)']);
    if k == 1, title('Velocity (ENU)'); legend(legendNames(hasVT),'Location','best'); end
end
subplot(4,1,4); hold on; grid on;
plot(tCT, sqrt(sum(navCT.usrVelENU.^2,2)), 'Color', cCT);
if hasVT, plot(tVT, sqrt(sum(navVT.usrVelENU.^2,2)), 'Color', cVT); end
ylabel('Speed (m/s)'); xlabel('Time since start of positioning (s)');
saveFig(f3, fullfile(outDir,'03_velocity_ENU.png'));

%% 4) C/N0
f4 = figure('Name','C/N0','Color','w','Position',[140 140 1000 750]);
nP = 1 + hasVT;
subplot(nP,1,1); plotCN0(CN0_CT, 0, svList); title('C/N0 - conventional tracking');
if hasVT
    subplot(nP,1,2); plotCN0(CN0_VT, file.skiptimeVT*1e-3, svList); title('C/N0 - vector tracking');
end
xlabel('Time since start of positioning (s)');
saveFig(f4, fullfile(outDir,'04_CN0.png'));

%% 5) Clock
f5 = figure('Name','Clock','Color','w','Position',[170 170 1000 500]);
dCT = navCT.clkDrift(:); 
subplot(2,1,1); hold on; grid on;
plot(tCT, navCT.clkBias(:) - navCT.clkBias(1), 'Color', cCT);
if hasVT, plot(tVT, navVT.clkBias(:) - navCT.clkBias(1), 'Color', cVT); end
ylabel('Clock bias change (m)'); title('Receiver clock'); legend(legendNames(hasVT),'Location','best');
subplot(2,1,2); hold on; grid on;
plot(tCT, dCT, 'Color', cCT);
if hasVT, plot(tVT, navVT.clkDrift(:), 'Color', cVT); end
ylabel('Clock drift (m/s)'); xlabel('Time since start of positioning (s)');
saveFig(f5, fullfile(outDir,'05_clock.png'));

%% Summary CSV (10 Hz)
step = 100;
T = [tCT(1:step:end) navCT.usrPosLLH(1:step:end,:) navCT.usrPosENU(1:step:end,:) navCT.usrVelENU(1:step:end,:)];
writeCsv(fullfile(outDir,'navCT_10Hz.csv'), T, 't_s,lat_deg,lon_deg,h_m,E_m,N_m,U_m,vE_mps,vN_mps,vU_mps');
if hasVT
    sV = max(1, round(0.1/median(diff(tVT))));            % 10 Hz
    T = [tVT(1:sV:end) navVT.usrPosLLH(1:sV:end,:) navVT.usrPosENU(1:sV:end,:) navVT.usrVelENU(1:sV:end,:)];
    writeCsv(fullfile(outDir,'navVT_10Hz.csv'), T, 't_s,lat_deg,lon_deg,h_m,E_m,N_m,U_m,vE_mps,vN_mps,vU_mps');
end
fprintf('Result figures and CSV files saved in %s\n', outDir);
end

function n = legendNames(hasVT)
if hasVT, n = {'CT (DLL/PLL + EKF)','VT (VDLL + EKF)'}; else, n = {'CT (DLL/PLL + EKF)'}; end
end

function plotCN0(CN0, t0, svList)
hold on; grid on;
nSv = size(CN0,2);
t   = t0 + (1:size(CN0,1))'*0.02;   % one C/N0 value per 20 ms
c   = lines(nSv);
names = cell(1,nSv);
for k = 1:nSv
    v = CN0(:,k); v(v==0) = NaN;
    plot(t, v, 'Color', c(k,:));
    names{k} = sprintf('PRN %d', svList(k));
end
ylabel('C/N0 (dB-Hz)'); ylim([20 60]);
legend(names,'Location','eastoutside');
end

function saveFig(f, name)
set(f,'PaperPositionMode','auto');
print(f, name, '-dpng', '-r110');
end

function writeCsv(name, M, header)
fid = fopen(name,'w'); fprintf(fid,'%s\n',header); fclose(fid);
dlmwrite(name, M, '-append', 'precision', '%.9g');
end
