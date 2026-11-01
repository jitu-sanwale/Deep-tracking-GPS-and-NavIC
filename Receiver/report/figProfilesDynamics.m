function figProfilesDynamics(sim, outDir)
%FIGPROFILESDYNAMICS Load factor and body roll rate of six maneuver profiles.
profDir = fullfile(sim, 'FighterProfiles');
names = {'turn','break_gonset','jinking','splitS_immelmann','rolls','departure_spin'};
ttl = {'Max-performance turn','Break turns','Jinking','Split-S / Immelmann','Rolls','Departure / spin'};
c1 = [0.165 0.471 0.839]; c2 = [0.922 0.408 0.204];
f = figure('Color','w','Units','inches','Position',[0.5 0.5 7.16 8.4],'Visible','off');
tl = tiledlayout(f, 6, 2, 'TileSpacing','tight','Padding','compact');
for k = 1:numel(names)
  S = load(fullfile(profDir, ['Fighter_' names{k} '_IMU_200Hz.mat']));
  F = S.imu.profile.flight;  w = S.imu.profile.maneuverWin;
  for j = 1:2
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    set(ax,'FontName','Times','FontSize',8,'GridAlpha',0.18,'TickDir','out', ...
        'LineWidth',1.1,'XColor',[0.1 0.1 0.1],'YColor',[0.1 0.1 0.1]);
    xregion(ax, w(1), w(2), 'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.12);
    if j == 1
      plot(ax, F.t, F.n, '-', 'Color', c1, 'LineWidth', 1.9); ylabel(ax,'n_z (g)');
      title(ax, [ttl{k} ': load factor'], 'FontWeight','normal');
    else
      plot(ax, F.t, F.wDeg(:,1), '-', 'Color', c2, 'LineWidth', 1.7); ylabel(ax,'p (deg/s)');
      title(ax, [ttl{k} ': roll rate'], 'FontWeight','normal');
    end
    xlim(ax,[0 60]);
    if k == numel(names), xlabel(ax,'Time (s)'); else, ax.XTickLabel = []; end
  end
end
exportgraphics(f, fullfile(outDir,'profiles_dynamics.png'), 'Resolution', 300);
close(f);
fprintf('written profiles_dynamics.png\n');
end
