function figAntennaAttitude(sim, outDir)
%FIGANTENNAATTITUDE Sketch of the body-fixed receive-antenna pattern
%   and the line of sight (LOS) to satellites, for level flight and a steep
%   bank, plus the gain versus body-frame elevation. Uses the simulator's
%   receive pattern (Simulator/Rx_GNSS_Antenna_Pattern.csv) and its
%   below-horizon roll-off (0.5 dB/deg, floor -20 dB), the attitude
%   transformation of applyAttitudeToBodyFrame.m, and the geometry of NavIC
%   PRN 2 and PRN 9 over Nashik at the start of the data.
%   Writes antenna_attitude.png (300 dpi) to outDir.

addpath(sim);
pat  = loadAntennaPattern(fullfile(sim, 'Rx_GNSS_Antenna_Pattern.csv'));
slope = 0.5; floorDb = -20;
gain = @(eb) arrayfun(@(e) lookupAntennaGain(pat, 0, e, slope, floorDb), eb);

heading = 90; pitch = 0;                     % flying east, level
rolls   = [0 85];                            % level, steep right bank
sats = struct('tag', {'PRN 2', 'PRN 9'}, 'name', {'NavIC PRN 2 (IGSO)', 'NavIC PRN 9 (IGSO)'}, ...
              'az', {300.5, 203.7}, 'el', {66.8, 30.8}, ...
              'col', {[0.80 0.15 0.12], [0.05 0.40 0.80]});
for s = 1:2
    for c = 1:2
        [~, sats(s).eb(c)] = applyAttitudeToBodyFrame(sats(s).az, sats(s).el, rolls(c), pitch, heading);
        sats(s).g(c) = gain(sats(s).eb(c));
    end
end
gmin = 3 + floorDb;

f = figure('Color','w','Units','centimeters','Position',[2 2 18 15.5],'Visible','off');
pos3d = {[0.01 0.43 0.47 0.55], [0.47 0.43 0.47 0.55]};
ttl   = {'(a) Level flight, heading east', sprintf('(b) %d%s right bank, heading east', rolls(2), char(176))};
for c = 1:2
    ax = axes(f, 'Position', pos3d{c});
    drawScene(ax, rolls(c), pitch, heading, sats, gain, gmin);
    title(ax, ttl{c}, 'FontWeight','normal', 'FontSize',10);
    % caption box with the LOS numbers
    lines = cell(1,2);
    for s = 1:2
        lines{s} = sprintf('\\color[rgb]{%.2f,%.2f,%.2f}%s:  e = %.1f%s,  e_b = %.1f%s,  G_{rx} = %+.1f dBi', ...
            sats(s).col, ['NavIC ' sats(s).tag], sats(s).el, char(176), sats(s).eb(c), char(176), sats(s).g(c));
    end
    annotation(f, 'textbox', [pos3d{c}(1)+0.02 0.415 pos3d{c}(3)-0.03 0.06], 'String', lines, ...
        'FontSize',7.5, 'FitBoxToText','off', 'EdgeColor',[0.25 0.25 0.25], 'LineWidth',1.2, 'BackgroundColor','w', 'Interpreter','tex', ...
        'HorizontalAlignment','left', 'VerticalAlignment','middle', 'Margin',3);
end
% shared colour bar for the gain lobe
cbAx = axes(f, 'Position',[0.955 0.52 0.001 0.38], 'Visible','off');
colormap(cbAx, turbo); clim(cbAx, [gmin 3]);
cb = colorbar(cbAx, 'Position',[0.945 0.52 0.014 0.38]);
cb.Label.String = 'G_{rx} (dBi)';  cb.FontSize = 8;

% (c) gain versus body-frame elevation
ax = axes(f, 'Position',[0.09 0.08 0.86 0.29]); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'FontSize',9,'TickDir','out','GridAlpha',0.15,'LineWidth',1.3,'XColor',[0.1 0.1 0.1],'YColor',[0.1 0.1 0.1]);
e = -90:0.25:90;  g = gain(e);
patch(ax, [-90 0 0 -90], [-24 -24 6 6], [0.94 0.94 0.94], 'EdgeColor','none');
text(ax, -45, 3.3, 'below the antenna horizon: back-lobe roll-off 0.5 dB/deg, floor -20 dB', ...
    'HorizontalAlignment','center','FontSize',7.5,'Color',[0.35 0.35 0.35]);
text(ax, 45, -12, 'above the antenna horizon: hemispherical pattern', ...
    'HorizontalAlignment','center','FontSize',7.5,'Color',[0.35 0.35 0.35]);
plot(ax, e, g, 'k-', 'LineWidth',2.6);
xline(ax, 0, ':', 'Color',[0.4 0.4 0.4]);
off = {[17 -1.2], [0 -3.6]};                 % label offsets [dx dy] for the banked points
for s = 1:2
    plot(ax, sats(s).eb, sats(s).g, ':', 'Color', sats(s).col, 'LineWidth',1.8);
    plot(ax, sats(s).eb(1), sats(s).g(1), 'o', 'MarkerSize',7, 'MarkerFaceColor','w', 'MarkerEdgeColor',sats(s).col, 'LineWidth',2.4);
    plot(ax, sats(s).eb(2), sats(s).g(2), 'o', 'MarkerSize',7, 'MarkerFaceColor',sats(s).col, 'MarkerEdgeColor',sats(s).col);
    text(ax, sats(s).eb(2)+off{s}(1), sats(s).g(2)+off{s}(2), sprintf('%s: %+.1f %s %+.1f dB', ...
        sats(s).tag, sats(s).g(1), char(8594), sats(s).g(2)), 'Color',sats(s).col, 'FontSize',8, ...
        'HorizontalAlignment','center');
end
xlim(ax, [-90 90]); ylim(ax, [-24 6]);
xlabel(ax, 'Body-frame elevation e_b of the LOS above the antenna horizon (deg)');
ylabel(ax, 'Receive gain G_{rx} (dBi)');
title(ax, '(c) Receive gain along the LOS  (open: level flight (a);  filled: bank (b))', 'FontWeight','normal','FontSize',10);

exportgraphics(f, fullfile(outDir, 'antenna_attitude.png'), 'Resolution', 300);
close(f);
fprintf('written %s\n', fullfile(outDir, 'antenna_attitude.png'));
end

%% ------------------------------------------------------------------------
function drawScene(ax, roll, pitch, heading, sats, gain, gmin)
hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
Cnb = dcmNED2Body(roll, pitch, heading);   % NED -> body
P   = @(xb) ned2plot(Cnb'*xb);             % body point -> plot (E, N, Up)

drawAircraft(ax, P);

% antenna on the spine, its horizon plane (body x-y plane) and gain lobe
ant = [0.12; 0; -0.07];
a0  = P(ant);
th = linspace(0, 2*pi, 97);
D = P(ant + 0.62*[cos(th); sin(th); zeros(size(th))]);
fill3(ax, D(1,:), D(2,:), D(3,:), [0.55 0.55 0.6], 'FaceAlpha',0.12, 'EdgeColor',[0.45 0.45 0.5], 'LineStyle','--', 'LineWidth',1.3);
[A, E] = meshgrid(linspace(0, 360, 97), linspace(-90, 90, 91));
G = reshape(gain(E(:)'), size(E));
r = 0.02 + 0.48*(G - gmin)/(3 - gmin);
X = zeros([size(E) 3]);
for i = 1:numel(E)
    [I, J] = ind2sub(size(E), i);
    X(I,J,:) = P(ant + r(i)*[cosd(E(i))*cosd(A(i)); cosd(E(i))*sind(A(i)); -sind(E(i))]);
end
surf(ax, X(:,:,1), X(:,:,2), X(:,:,3), G, 'EdgeColor','none', 'FaceAlpha',0.45, ...
    'FaceLighting','gouraud', 'AmbientStrength',0.6);
colormap(ax, turbo); clim(ax, [gmin 3]);
bs = P(ant + [0;0;-0.72]);
plot3(ax, [a0(1) bs(1)], [a0(2) bs(2)], [a0(3) bs(3)], 'k--', 'LineWidth',1.4);
text(ax, bs(1), bs(2), bs(3), '  boresight -z_b', 'FontSize',7.5);

% body axes from the c.g.
L = 0.55; axd = {[L;0;0], 'x_b', [0.70 0.05 0.05]; [0;L;0], 'y_b', [0.05 0.50 0.10]; [0;0;L], 'z_b', [0.10 0.20 0.75]};
o = P([0;0;0]);
for k = 1:3
    t = P(axd{k,1});
    quiver3(ax, o(1),o(2),o(3), t(1)-o(1),t(2)-o(2),t(3)-o(3), 0, 'Color',axd{k,3}, 'LineWidth',2.9, 'MaxHeadSize',0.45);
    q = o + 1.15*(t - o);
    text(ax, q(1), q(2), q(3), axd{k,2}, 'Color',axd{k,3}, 'FontSize',9, 'FontWeight','bold', 'HorizontalAlignment','center');
end

% LOS to the satellites, with the body-frame elevation arc
for s = 1:numel(sats)
    uN = [cosd(sats(s).el)*cosd(sats(s).az); cosd(sats(s).el)*sind(sats(s).az); -sind(sats(s).el)];
    p1 = a0 + 1.22*ned2plot(uN);
    quiver3(ax, a0(1),a0(2),a0(3), p1(1)-a0(1),p1(2)-a0(2),p1(3)-a0(3), 0, 'Color',sats(s).col, 'LineWidth',3.2, 'MaxHeadSize',0.2);
    plot3(ax, p1(1), p1(2), p1(3), 'p', 'MarkerSize',11, 'MarkerFaceColor',sats(s).col, 'MarkerEdgeColor','k');
    text(ax, p1(1), p1(2), p1(3)+0.09, sats(s).tag, 'Color',sats(s).col, 'FontSize',8.5, 'FontWeight','bold', 'HorizontalAlignment','center');
    % arc from the LOS projection on the antenna horizon plane up/down to the LOS
    ub = Cnb*uN;  hb = [ub(1); ub(2); 0]/norm(ub(1:2));
    eb = asin(-ub(3));
    tt = linspace(0, eb, 30);
    arc = ant + 0.34*(hb*cos(tt) + [0;0;-1]*sin(tt));
    Q = P(arc);  H = P(ant + 0.62*hb);
    plot3(ax, [a0(1) H(1)], [a0(2) H(2)], [a0(3) H(3)], ':', 'Color',sats(s).col, 'LineWidth',1.6);
    plot3(ax, Q(1,:), Q(2,:), Q(3,:), '-', 'Color',sats(s).col, 'LineWidth',2.1);
    m = Q(:, round(end/2));
    text(ax, m(1), m(2), m(3), ' e_b', 'Color',sats(s).col, 'FontSize',8);
end
view(ax, [-40 16]);
camlight(ax, 'headlight');
lim = 1.15; xlim(ax, [-lim lim]); ylim(ax, [-lim lim]); zlim(ax, [-0.75 1.35]);
end

function drawAircraft(ax, P)
% simple fighter-like shape in body axes (x forward, y right, z down)
grey = [0.68 0.70 0.74];
xs = linspace(-0.55, 0.58, 50);
rr = 0.06*sqrt(max(1 - (xs/0.59).^2, 0)).*(1 - 0.3*(xs < -0.32));
ph = linspace(0, 2*pi, 33);
[XX, PH] = meshgrid(xs, ph);  RR = repmat(rr, numel(ph), 1);
Xf = zeros([size(XX) 3]);
for i = 1:numel(XX)
    [I, J] = ind2sub(size(XX), i);
    Xf(I,J,:) = P([XX(i); RR(i)*cos(PH(i)); RR(i)*sin(PH(i))]);
end
surf(ax, Xf(:,:,1), Xf(:,:,2), Xf(:,:,3), 'FaceColor',grey, 'EdgeColor','none', 'FaceLighting','gouraud');
for sgn = [1 -1]
    fillB(ax, P, [0.20 sgn*0.05 0; -0.17 sgn*0.50 0; -0.28 sgn*0.50 0; -0.26 sgn*0.05 0]', grey);   % wing
    fillB(ax, P, [-0.36 sgn*0.04 0; -0.49 sgn*0.24 0; -0.55 sgn*0.24 0; -0.55 sgn*0.04 0]', grey);  % stabilator
end
fillB(ax, P, [-0.28 0 -0.055; -0.50 0 -0.32; -0.56 0 -0.32; -0.56 0 -0.055]', grey);              % fin
a = P([0.12; 0; -0.062]);
plot3(ax, a(1), a(2), a(3), 's', 'MarkerSize',5, 'MarkerFaceColor',[0.95 0.75 0.1], 'MarkerEdgeColor','k');
end

function fillB(ax, P, V, col)
W = P(V);
fill3(ax, W(1,:), W(2,:), W(3,:), col, 'EdgeColor',[0.3 0.3 0.3], 'LineWidth',1.0, 'FaceLighting','gouraud');
end

function q = ned2plot(n)
q = [n(2,:); n(1,:); -n(3,:)];   % NED -> (E, N, Up)
end

function C = dcmNED2Body(roll, pitch, heading)
sp = sind(roll); cp = cosd(roll); st = sind(pitch); ct = cosd(pitch); ss = sind(heading); cs = cosd(heading);
C = [ct*cs, ct*ss, -st; -cp*ss + sp*st*cs, cp*cs + sp*st*ss, sp*ct; sp*ss + cp*st*cs, -sp*cs + cp*st*ss, cp*ct];
end
