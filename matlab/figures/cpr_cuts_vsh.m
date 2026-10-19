function info = cpr_cuts_vsh(cutdir, ratios, k)
%CPR_CUTS_VSH  Truth vs. recovered pattern cut through the pattern peak, with the pointwise error, for the
%phaseless IS examples of pr_examples with the spherical-wave basis and both polarizations (paper Fig. 5,
%supplementary azimuth cut). Dual-polarized version of cpr_cuts.m (2026-10-19).
%   cpr_cuts_vsh()                 elevation cut, examples at m/n = 0.10 and 0.30, first frequency
%   cpr_cuts_vsh('phi')            azimuth cut through the peak
% The truth is the exact pattern (pattern_crlh_vsh). The plotted component is the one that holds the peak
% of the pattern at frequency k (along or across the scan circle). One global phase over all wavelengths and
% both polarizations is removed (the scoring convention for K > 1). Four stacked panels: magnitude, magnitude
% of the complex error relative to the peak ("|Error| (dB)"), phase, phase error ("Error (deg)").
% Writes <results>/figs/paper/pr3_cut.pdf or pr3_cut_phi.pdf.
if nargin < 1 || isempty(cutdir), cutdir = 'theta'; end
if nargin < 2 || isempty(ratios), ratios = [0.10 0.30]; end
if nargin < 3 || isempty(k), k = 1; end
outDir = cpr_paths().results; pDir = fullfile(outDir, 'figs', 'paper');
if ~exist(pDir, 'dir'), mkdir(pDir); end
set(0,'defaultAxesFontName','Times'); set(0,'defaultTextFontName','Times');
set(0,'defaultTextFontSize',11);      set(0,'defaultAxesFontSize',11); set(0,'defaultLegendFontName','Times');
E = {};
for d = dir(fullfile(outDir, 'pr_examples_w*.mat')).'
    S0 = load(fullfile(d.folder, d.name));
    for e = S0.res(:).', if ~isempty(e.Xh), E{end+1} = e; end, end %#ok<AGROW>
end
rt = cellfun(@(e) e.cfg.ratio, E);
sel = arrayfun(@(r) find(abs(rt - r) < 1e-9, 1), ratios); E = E(sel);
c = E{1}.cfg; [X, ~] = pattern_crlh_vsh(c.nw, c.nd, c.f0idx, 16, c.npol);
nd = c.nd; nh = c.nh; nw = c.nw;
for i = 1:numel(E)
    Xh = E{i}.Xh; E{i}.Xh = Xh * exp(1i*angle(Xh(:)'*X(:)));
end
% component holding the peak at frequency k
pk = arrayfun(@(q) max(max(abs(X(:, :, (q-1)*nd + k)))), 1:c.npol); [~, q] = max(pk);
s = (q-1)*nd + k; US = X(:, :, s); peakMag = max(abs(US(:)));
[~, ip] = max(abs(US(:))); [ith, iph] = ind2sub(size(US), ip);
errdB = @(a, b) max(20*log10(abs(a - b) / peakMag), -100);
wrapdeg = @(a) angle(exp(1i*a)) * 180/pi;
if strcmp(cutdir, 'phi')
    ith2 = mod(nw - (ith - 1), nw) + 1;
    if (ith - 1)*360/nw > 180, [ith, ith2] = deal(ith2, ith); end
    cut = @(M) [M(ith, :), M(ith2, :)].';
    C.theta = (0:2*nh-1)*180/nh; C.thetaLim = [0 360-180/nh]; C.xlab = '$\phi$ ($^\circ$)';
    outName = 'pr3_cut_phi'; info.cutAt = sprintf('polar angle %.2f deg', (ith-1)*360/nw);
else
    cut = @(M) M(:, iph);
    C.theta = (0:nw-1)*360/nw; C.thetaLim = [0 360-360/nw]; C.xlab = '$\theta$ ($^\circ$)';
    outName = 'pr3_cut'; info.cutAt = sprintf('phi %.2f deg', (iph-1)*180/nh);
end
elT = cut(US); elR = cellfun(@(e) cut(double(e.Xh(:, :, s))), E, 'UniformOutput', false);
mn = char(8722);
C.nm = cellfun(@(e) sprintf('%s%.1f dB', mn, -e.global), E, 'UniformOutput', false);
C.mT = 20*log10(abs(elT)/peakMag);
C.mR = cellfun(@(v) 20*log10(abs(v)/peakMag), elR, 'UniformOutput', false);
C.mE = cellfun(@(v) errdB(v, elT), elR, 'UniformOutput', false);
C.pT = angle(elT)*180/pi;
C.pR = cellfun(@(v) angle(v)*180/pi, elR, 'UniformOutput', false);
C.pE = cellfun(@(v) wrapdeg(angle(v) - angle(elT)), elR, 'UniformOutput', false);
make_fig4(C, fullfile(pDir, outName));
comp = {'along the scan circle', 'across the scan circle'};
info.component = comp{q}; info.nmse = cellfun(@(e) e.global, E); info.ratios = ratios;
info.freqNMSE = cellfun(@(e) 20*log10(norm(reshape(e.Xh(:,:,[k, nd+k]) - X(:,:,[k, nd+k]), [], 1)) / ...
    norm(reshape(X(:,:,[k, nd+k]), [], 1))), E);
fprintf('%s: component %s, cut at %s, NMSE %s dB (whole), %s dB at frequency %d\n', outName, info.component, ...
    info.cutAt, mat2str(round(info.nmse, 1)), mat2str(round(info.freqNMSE, 1)), k);
end

function make_fig4(S, outBase)
W = 12; H = 10.5; left = 1.6; right = 0.1; top = 0.2; bottom = 1.0; gap = 0.3; legGap = 0.25;
lw = 1; lwRec = 1.6;
cTruth = [0 0.4470 0.7410]; cRec = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250]; ls = {'--', '-.'};
fig = figure('Visible','on','Units','centimeters','Position',[2 2 W H], 'Color','w');
set(fig, 'PaperUnits','centimeters', 'PaperSize',[W H], 'PaperPosition',[0 0 W H]);
boxH = (H - top - bottom - 3*gap) / 4;
y = @(k) bottom + (4 - k)*(boxH + gap);
w0 = W - left - right - 3;
ax = gobjects(1, 4);
for k = 1:4, ax(k) = axes('Units','centimeters','Position',[left, y(k), w0, boxH]); hold(ax(k), 'on'); end
hT = plot(ax(1), S.theta, S.mT, 'Linewidth', lw, 'color', cTruth);
hP = plot(ax(3), S.theta, S.pT, 'Linewidth', lw, 'color', cTruth);
hR = gobjects(1, numel(S.mR));
for i = 1:numel(S.mR)
    hR(i) = plot(ax(1), S.theta, S.mR{i}, ls{i}, 'Linewidth', lwRec, 'color', cRec(i, :));
    plot(ax(2), S.theta, S.mE{i}, ls{i}, 'Linewidth', 1.2, 'color', cRec(i, :));
    plot(ax(3), S.theta, S.pR{i}, ls{i}, 'Linewidth', lwRec, 'color', cRec(i, :));
    plot(ax(4), S.theta, S.pE{i}, ls{i}, 'Linewidth', 1.2, 'color', cRec(i, :));
end
uistack(hT, 'top'); uistack(hP, 'top');
ylabel(ax(1), '$\left| E \right|$ (dB)', 'Interpreter', 'latex');
ylabel(ax(2), '$\left|\mathrm{Error}\right|$ (dB)', 'Interpreter', 'latex');
ylabel(ax(3), '$\angle E$ ($^\circ$)', 'Interpreter', 'latex');
ylabel(ax(4), 'Error ($^\circ$)', 'Interpreter', 'latex');
ylim(ax(3), [-180 180]); ax(3).YTick = -180:90:180;
for k = 1:4, xlim(ax(k), S.thetaLim); box(ax(k), 'on'); end
for k = 1:3, ax(k).XTickLabel = []; end
for k = [1 2 4], more_yticks(ax(k)); end
xlabel(ax(4), S.xlab, 'Interpreter', 'latex');
lg = legend(ax(1), [hT hR], [{'Truth'}, S.nm]); drawnow
lg.Units = 'centimeters'; legW = 1.08*lg.Position(3);
boxW = W - left - right - legW - legGap;
for k = 1:4, ax(k).Position = [left, y(k), boxW, boxH]; end
lg.Position(3) = legW; lg.Position(1:2) = [left + boxW + legGap, y(1) + boxH - lg.Position(4)];
saveas(fig, outBase, 'pdf'); fprintf('wrote %s.pdf\n', outBase);
end

function more_yticks(ax)
drawnow; yl = ax.YLim;
for s = [1 2 5 10 20 50 100 0.5 4]
    t = ceil(yl(1)/s - 1e-9)*s : s : yl(2) + 1e-9;
    if numel(t) >= 4 && numel(t) <= 6, ax.YTick = t; return; end
    lo = floor(yl(1)/s)*s; hi = ceil(yl(2)/s)*s;
    n = round((hi - lo)/s) + 1;
    if n >= 4 && n <= 6 && yl(1) - lo < s/3 && hi - yl(2) < s/3
        ax.YLim = [lo hi]; ax.YTick = lo:s:hi; return;
    end
end
end
