function cpr_cut_panels(S, outBase, showTruthPhaseOnTop)
%CPR_CUT_PANELS  2x2 cut figure (2026-10-25): pattern cuts on top (magnitude left, phase right), the pointwise errors below
%(|error| re the peak left, phase error right), one legend row above. Used by cpr_cuts_vsh and cpr_fwc_cuts.
%   S: theta, thetaLim, xlab, nm (cell of labels), mT, mR{}, mE{}, pT, pR{}, pE{}  (as the four-panel version)
% Size 9 x 7.2 cm, so two of them fit side by side in a two-column figure. Writes <outBase>.pdf and .png.
if nargin < 3, showTruthPhaseOnTop = true; end %#ok<NASGU>
W = 9; H = 7.0; left = 1.5; mid = 1.6; right = 0.15; top = 0.9; bottom = 1.2; gap = 0.35; fs = 10;
lw = 1; lwRec = 1.4;
cTruth = [0 0.4470 0.7410]; cRec = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250]; ls = {'--', '-.'};
fig = figure('Visible','off','Units','centimeters','Position',[2 2 W H], 'Color','w');
set(fig, 'PaperUnits','centimeters', 'PaperSize',[W H], 'PaperPosition',[0 0 W H]);
bw = (W - left - mid - right)/2; bh = (H - top - bottom - gap)/2;
x = @(c) left + (c-1)*(bw + mid); y = @(r) bottom + (2 - r)*(bh + gap);
ax = gobjects(2, 2);
for r = 1:2, for c = 1:2
    ax(r, c) = axes('Units','centimeters','Position',[x(c), y(r), bw, bh], 'FontName','Times','FontSize',fs); hold(ax(r, c), 'on');
end, end
hT = plot(ax(1, 1), S.theta, S.mT, 'Linewidth', lw, 'color', cTruth);
hP = plot(ax(1, 2), S.theta, S.pT, 'Linewidth', lw, 'color', cTruth);
hR = gobjects(1, numel(S.mR));
for i = 1:numel(S.mR)
    hR(i) = plot(ax(1, 1), S.theta, S.mR{i}, ls{i}, 'Linewidth', lwRec, 'color', cRec(i, :));
    plot(ax(2, 1), S.theta, S.mE{i}, ls{i}, 'Linewidth', 1.1, 'color', cRec(i, :));
    plot(ax(1, 2), S.theta, S.pR{i}, ls{i}, 'Linewidth', lwRec, 'color', cRec(i, :));
    plot(ax(2, 2), S.theta, S.pE{i}, ls{i}, 'Linewidth', 1.1, 'color', cRec(i, :));
end
uistack(hT, 'top'); uistack(hP, 'top');
ylabel(ax(1, 1), '$\left| E \right|$ (dB)', 'Interpreter', 'latex', 'FontSize', fs);
ylabel(ax(2, 1), '$\left|\mathrm{Error}\right|$ (dB)', 'Interpreter', 'latex', 'FontSize', fs);
ylabel(ax(1, 2), '$\angle E$ ($^\circ$)', 'Interpreter', 'latex', 'FontSize', fs);
ylabel(ax(2, 2), 'Error ($^\circ$)', 'Interpreter', 'latex', 'FontSize', fs);
ylim(ax(1, 2), [-180 180]); ax(1, 2).YTick = -180:180:180;
for r = 1:2, for c = 1:2, xlim(ax(r, c), S.thetaLim); box(ax(r, c), 'on'); end, end
for c = 1:2
    ax(1, c).XTickLabel = [];
    xlabel(ax(2, c), S.xlab, 'Interpreter', 'latex', 'FontSize', fs);
end
for a = ax(:).', a.XTick = 0:180:360; a.XTickLabelRotation = 0; end
more_yticks(ax(1, 1)); more_yticks(ax(2, 1)); more_yticks(ax(2, 2));
lg = legend(ax(1, 1), [hT hR], [{'Truth'}, S.nm], 'Orientation', 'horizontal', 'FontName', 'Times', 'FontSize', fs); drawnow
lg.Units = 'centimeters'; lg.Position(1:2) = [(W - lg.Position(3))/2, H - lg.Position(4) - 0.08];
saveas(fig, outBase, 'pdf'); print(fig, outBase, '-dpng', '-r150'); fprintf('wrote %s.pdf\n', outBase);
close(fig);
end

function more_yticks(ax)
drawnow; yl = ax.YLim;
for s = [1 2 5 10 20 50 100 0.5 4]
    t = ceil(yl(1)/s - 1e-9)*s : s : yl(2) + 1e-9;
    if numel(t) >= 3 && numel(t) <= 5, ax.YTick = t; return; end
    lo = floor(yl(1)/s)*s; hi = ceil(yl(2)/s)*s;
    n = round((hi - lo)/s) + 1;
    if n >= 3 && n <= 5 && yl(1) - lo < s/3 && hi - yl(2) < s/3
        ax.YLim = [lo hi]; ax.YTick = lo:s:hi; return;
    end
end
end
