function info = cpr_fwc_cuts(fGHz, thetaCut, phiCut)
%CPR_FWC_CUTS  Section VII cut figures from the full-ring trials (fwc_run_trials -> fwc_trials.mat, 2026-10-21): free-space
%far field of the AUT, truth vs. the recovered coefficients of the WORST of the trials, with the pointwise error, in the
%four-panel format of cpr_cuts_vsh (magnitude, |error| re peak, phase, phase error).
%   cpr_fwc_cuts()                     2.725 GHz, azimuth cut at theta = 90 deg, elevation cut at phi = 0 deg
%   cpr_fwc_cuts(fGHz, thetaCut, phiCut)
% The recovered coefficients are aligned with one complex factor at that frequency (the main score of Section VII).
% Far field: c = a .* j^(l+1) (TE) or j^l (TM) (fwc_aut_modes), E = sum c (A1 or A2) (vsh_eval). The plotted component is
% the one with more energy along each cut (theta-hat or phi-hat). The elevation cut runs over the full great circle
% through both poles (theta from 0 to 360 deg; beyond 180 deg it is the opposite half-plane, theta-hat reversed), so it
% includes the polar caps that no probe measures.
% Reads <results>/fwc_validation.mat and fwc_trials.mat; writes <results>/figs/paper/cutphi.pdf and cuttheta.pdf.
if nargin < 1 || isempty(fGHz), fGHz = 2.725; end
if nargin < 2 || isempty(thetaCut), thetaCut = 90; end
if nargin < 3 || isempty(phiCut), phiCut = 0; end
outDir = cpr_paths().results; pDir = fullfile(outDir, 'figs', 'paper');
if ~exist(pDir, 'dir'), mkdir(pDir); end
vFile = fullfile(outDir, 'fwc_validation.mat');
if ~exist(vFile, 'file'), vFile = fullfile(cpr_paths().data, 'fullwave', 'fwc_validation.mat'); end   % released data
S = load(vFile); V = S.V;
Tr = load(fullfile(outDir, 'fwc_trials.mat')); T = Tr.T;
sc = arrayfun(@(t) t.nmse.perfreqScale, T); [~, iw] = max(sc); tr = T(iw);
[~, k] = min(abs(V.freq - fGHz*1e9));
modes = V.A.modes; L = max(modes(:, 1));
fac = 1i .^ (modes(:, 1) + (modes(:, 3) == 1));
cT = V.A.a(:, k) .* fac; cR = tr.aHat(:, k) * tr.perfreqFactor(k) .* fac;
set(0,'defaultAxesFontName','Times'); set(0,'defaultTextFontName','Times');
set(0,'defaultTextFontSize',11);      set(0,'defaultAxesFontSize',11); set(0,'defaultLegendFontName','Times');
mn = char(8722); nm = {sprintf('%s%.1f dB', mn, -tr.nmse.perfreqScale)};
% azimuth cut
ph = (0:1:359).' * pi/180; th = thetaCut*pi/180 * ones(size(ph));
[Yt, Yp] = vsh_eval(L, th, ph);
[uT, uR, cn] = pick(Yt*cT, Yp*cT, Yt*cR, Yp*cR);
C = panels(uT, uR, nm); C.theta = ph*180/pi; C.thetaLim = [0 359]; C.xlab = '$\phi$ ($^\circ$)';
make_fig4(C, fullfile(pDir, 'cutphi'));
info.azimuth = sprintf('theta = %g deg, component %s', thetaCut, cn);
% elevation cut over the full great circle
s = (0.5:1:359.5).'; far = s > 180; th = s; th(far) = 360 - s(far); pp = phiCut*ones(size(s)); pp(far) = phiCut + 180;
[Yt, Yp] = vsh_eval(L, th*pi/180, pp*pi/180); sg = ones(size(s)); sg(far) = -1;
[uT, uR, cn] = pick(sg.*(Yt*cT), Yp*cT, sg.*(Yt*cR), Yp*cR);
C = panels(uT, uR, nm); C.theta = s; C.thetaLim = [0 360]; C.xlab = '$\theta$ ($^\circ$)';
make_fig4(C, fullfile(pDir, 'cuttheta'));
info.elevation = sprintf('phi = %g deg (great circle), component %s', phiCut, cn);
info.seed = tr.seed; info.nmse = tr.nmse.perfreqScale; info.f = V.freq(k);
info.freqNMSE = 20*log10(norm(cR./fac - V.A.a(:, k)) / norm(V.A.a(:, k)));
fprintf('worst trial: seed %d, per-frequency score %.1f dB; f = %.3f GHz, coefficient NMSE there %.1f dB\n', ...
    info.seed, info.nmse, info.f/1e9, info.freqNMSE);
fprintf('cutphi: %s; cuttheta: %s\n', info.azimuth, info.elevation);
end

function [uT, uR, cn] = pick(tT, pT, tR, pR)
if norm(tT) >= norm(pT), uT = tT; uR = tR; cn = 'E_theta'; else, uT = pT; uR = pR; cn = 'E_phi'; end
end

function C = panels(uT, uR, nm)
peak = max(abs(uT));
C.nm = nm;
C.mT = 20*log10(abs(uT)/peak); C.mR = {20*log10(abs(uR)/peak)};
C.mE = {max(20*log10(abs(uR - uT)/peak), -100)};
C.pT = angle(uT)*180/pi; C.pR = {angle(uR)*180/pi};
C.pE = {angle(exp(1i*(angle(uR) - angle(uT))))*180/pi};
end

function make_fig4(S, outBase)
% four stacked panels, as cpr_cuts_vsh.m
W = 12; H = 10.5; left = 1.6; right = 0.1; top = 0.2; bottom = 1.0; gap = 0.3; legGap = 0.25;
lw = 1; lwRec = 1.6;
cTruth = [0 0.4470 0.7410]; cRec = [0.8500 0.3250 0.0980]; ls = {'--'};
fig = figure('Visible','off','Units','centimeters','Position',[2 2 W H], 'Color','w');
set(fig, 'PaperUnits','centimeters', 'PaperSize',[W H], 'PaperPosition',[0 0 W H]);
boxH = (H - top - bottom - 3*gap) / 4;
y = @(k) bottom + (4 - k)*(boxH + gap);
w0 = W - left - right - 3;
ax = gobjects(1, 4);
for k = 1:4, ax(k) = axes('Units','centimeters','Position',[left, y(k), w0, boxH]); hold(ax(k), 'on'); end
hT = plot(ax(1), S.theta, S.mT, 'Linewidth', lw, 'color', cTruth);
hP = plot(ax(3), S.theta, S.pT, 'Linewidth', lw, 'color', cTruth);
hR = plot(ax(1), S.theta, S.mR{1}, ls{1}, 'Linewidth', lwRec, 'color', cRec);
plot(ax(2), S.theta, S.mE{1}, ls{1}, 'Linewidth', 1.2, 'color', cRec);
plot(ax(3), S.theta, S.pR{1}, ls{1}, 'Linewidth', lwRec, 'color', cRec);
plot(ax(4), S.theta, S.pE{1}, ls{1}, 'Linewidth', 1.2, 'color', cRec);
uistack(hT, 'top'); uistack(hP, 'top');
ylabel(ax(1), '$\left| E \right|$ (dB)', 'Interpreter', 'latex');
ylabel(ax(2), '$\left|\mathrm{Error}\right|$ (dB)', 'Interpreter', 'latex');
ylabel(ax(3), '$\angle E$ ($^\circ$)', 'Interpreter', 'latex');
ylabel(ax(4), 'Error ($^\circ$)', 'Interpreter', 'latex');
ylim(ax(3), [-180 180]); ax(3).YTick = -180:90:180;
for k = 1:4, xlim(ax(k), S.thetaLim); box(ax(k), 'on'); end
for k = 1:3, ax(k).XTickLabel = []; end
xlabel(ax(4), S.xlab, 'Interpreter', 'latex');
lg = legend(ax(1), [hT hR], [{'Truth'}, S.nm]); drawnow
lg.Units = 'centimeters'; legW = 1.08*lg.Position(3);
boxW = W - left - right - legW - legGap;
for k = 1:4, ax(k).Position = [left, y(k), boxW, boxH]; end
lg.Position(3) = legW; lg.Position(1:2) = [left + boxW + legGap, y(1) + boxH - lg.Position(4)];
saveas(fig, outBase, 'pdf'); print(fig, outBase, '-dpng', '-r150'); fprintf('wrote %s.pdf\n', outBase);
close(fig);
end
