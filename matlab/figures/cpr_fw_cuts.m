function cpr_fw_cuts(tag, seed, fGHz)
%CPR_FW_CUTS  Section VI cut figures from the full-wave trials: truth vs. recovered azimuth cut
%(theta = 57.3 deg, cutphi.pdf) and elevation cut (phi = 90 deg, cuttheta.pdf), each with the
%pointwise error below (12 cm wide, Times 11 pt).
%   cpr_fw_cuts('calop', 9, 2.725)   trial with seed 9 at the frequency nearest 2.725 GHz
%One global phase is removed over the whole pattern (K = 4 readings measure the phases between
%frequencies). Writes <results>/figs/paper/cutphi.pdf, cuttheta.pdf and cutEfull3/cutHfull3.pdf.
outDir = cpr_paths().results;
pDir = fullfile(outDir, 'figs', 'paper');
L = load(fullfile(outDir, 'fw_data_hfss.mat')); D = L.D;
S = load(fullfile(outDir, sprintf('fw_trials_%s.mat', tag))); R = S.R;
r = R([R.seed] == seed);
nf = size(r.Xh, 3);                                        % frequencies used by the trial (quick preset: fewer)
D.Xring = D.Xring(:, :, 1:nf); D.f = D.f(1:nf);
[~, t] = min(abs(D.f - fGHz*1e9));
X = D.Xring; Xh = double(r.Xh);
Xh = Xh * exp(1i*angle(Xh(:)'*X(:)));                    % one global phase
fprintf('%s seed %d: global NMSE %.2f dB; plotted f = %.3f GHz, NMSE there %.2f dB\n', tag, seed, ...
    r.score.global, D.f(t)/1e9, 20*log10(norm(reshape(Xh(:,:,t)-X(:,:,t),[],1))/norm(reshape(X(:,:,t),[],1))));
US = X(:,:,t); X1 = Xh(:,:,t);
[nw, nh] = size(US);
theta_deg = linspace(20, 90, nw); phi_deg = (0:nh-1)*360/nh;
slicep = round(nw/2) + 1;   slicep1 = round(nh/4) + 1;   % theta = 57.3 deg, phi = 90 deg
fprintf('azimuth cut theta = %.2f deg, elevation cut phi = %.2f deg\n', theta_deg(slicep), phi_deg(slicep1));
peakMag = max(abs(US(:)));
errdB = @(a, b) max(20*log10(abs(a - b) / peakMag), -100);
wrapdeg = @(a) angle(exp(1i*a)) * 180/pi;
el_truth = US(:,slicep1); el_rec = X1(:,slicep1);
az_truth = US(slicep,:).'; az_rec = X1(slicep,:).';
set(0,'defaultAxesFontName','Times'); set(0,'defaultTextFontName','Times');
set(0,'defaultTextFontSize',11); set(0,'defaultAxesFontSize',11);
P = struct('theta', theta_deg, 'phi', phi_deg, 'thetaLim', [theta_deg(1) theta_deg(end)], 'phiLim', [0 360-360/nh]);
P.elT = 20*log10(abs(el_truth)); P.elR = 20*log10(abs(el_rec));
P.azT = 20*log10(abs(az_truth)); P.azR = 20*log10(abs(az_rec));
P.elE = errdB(el_rec, el_truth); P.azE = errdB(az_rec, az_truth);
P.ylab1 = '$\left| \widetilde{E} \right|$ (dB)'; P.ylab2 = '$\left| \mathrm{Error} \right|$ (dB)';
make_fig(P, fullfile(pDir, 'cutEfull3'));
fprintf('max magnitude error along the cuts %.1f dB re peak\n', max([P.elE; P.azE]));
P.elT = angle(el_truth)*180/pi; P.elR = angle(el_rec)*180/pi;
P.azT = angle(az_truth)*180/pi; P.azR = angle(az_rec)*180/pi;
P.elE = wrapdeg(angle(el_rec) - angle(el_truth)); P.azE = wrapdeg(angle(az_rec) - angle(az_truth));
P.ylab1 = '$\angle E$ ($^\circ$)'; P.ylab2 = 'Phase error ($^\circ$)';
make_fig(P, fullfile(pDir, 'cutHfull3'));
fprintf('max |phase error| along the cuts %.1f deg\n', max(abs([P.elE; P.azE])));
% --- phi cut only: magnitude and phase vs. phi, two stacked panels
eM = errdB(az_rec, az_truth); eP = wrapdeg(angle(az_rec) - angle(az_truth));
fprintf('phi cut (theta = %.1f deg): max magnitude error %.1f dB re peak, max |phase error| %.1f deg\n', ...
    theta_deg(slicep), max(eM), max(abs(eP)));
make_phi_fig(phi_deg, az_truth, az_rec, peakMag, fullfile(pDir, 'cutphi'));
% --- theta cut (supplement): elevation cut across the probe arc at phi = 90 deg
make_phi_fig(theta_deg, el_truth, el_rec, peakMag, fullfile(pDir, 'cuttheta'), '$\theta$ ($^\circ$)', [theta_deg(1) theta_deg(end)]);
fprintf('theta cut (phi = %.2f deg): max magnitude error %.1f dB re peak, max |phase error| %.1f deg\n', ...
    phi_deg(slicep1), max(errdB(el_rec, el_truth)), max(abs(wrapdeg(angle(el_rec) - angle(el_truth)))));
end

function make_phi_fig(phi, xt, xr, pk, outBase, xlab, xl)
% phi: the cut's abscissa (phi by default; theta for the supplement's elevation cut)
if nargin < 6, xlab = '$\phi$ ($^\circ$)'; xl = [0 360 - 360/numel(phi)]; end
% four stacked panels: magnitude, magnitude error, phase, phase error
W = 12; H = 9.5; left = 1.6; right = 0.15; top = 0.2; bottom = 1.0; gap = 0.3;
boxH = (H - top - bottom - 3*gap) / 4; boxW = W - left - right;
lw = 1; lwRec = 1.6; cT = [0 0.4470 0.7410]; cR = [0.9290 0.6940 0.1250]; cE = [0.8500 0.3250 0.0980];
fig = figure('Units','centimeters','Position',[2 2 W H], 'Color','w', 'Visible', 'off');
set(fig, 'PaperUnits','centimeters', 'PaperSize',[W H], 'PaperPosition',[0 0 W H]);
y = @(k) bottom + (4 - k)*(boxH + gap);
a1 = axes('Units','centimeters','Position',[left, y(1), boxW, boxH]); hold(a1, 'on');
hR1 = plot(a1, phi, 20*log10(abs(xr)), '--', 'Linewidth', lwRec, 'color', cR);
hT1 = plot(a1, phi, 20*log10(abs(xt)), 'Linewidth', lw, 'color', cT);
ylabel(a1, '$|E|$ (dB)', 'Interpreter','latex');
legend(a1, [hT1 hR1], {'Truth','Recovered'}, 'Location', 'best');
a2 = axes('Units','centimeters','Position',[left, y(2), boxW, boxH]); hold(a2, 'on');
plot(a2, phi, max(20*log10(abs(xr - xt)/pk), -100), 'Linewidth', lw, 'color', cE);
ylabel(a2, '$|\mathrm{Error}|$ (dB)', 'Interpreter','latex');
a3 = axes('Units','centimeters','Position',[left, y(3), boxW, boxH]); hold(a3, 'on');
plot(a3, phi, angle(xr)*180/pi, '--', 'Linewidth', lwRec, 'color', cR);
plot(a3, phi, angle(xt)*180/pi, 'Linewidth', lw, 'color', cT);
ylabel(a3, '$\angle E$ ($^\circ$)', 'Interpreter','latex'); ylim(a3, [-180 180]); a3.YTick = -180:90:180;
a4 = axes('Units','centimeters','Position',[left, y(4), boxW, boxH]); hold(a4, 'on');
plot(a4, phi, angle(exp(1i*(angle(xr) - angle(xt))))*180/pi, 'Linewidth', lw, 'color', cE);
ylabel(a4, 'Error ($^\circ$)', 'Interpreter','latex');   % short: fits the shorter panels
xlabel(a4, xlab, 'Interpreter','latex');
for ax = [a1 a2 a3 a4], xlim(ax, xl); grid(ax, 'on'); box(ax, 'on'); end
for ax = [a1 a2 a3], ax.XTickLabel = []; end
for ax = [a1 a2 a4], more_yticks(ax); end            % 4-6 labels per panel
saveas(fig, outBase, 'pdf'); close(fig);
fprintf('wrote %s.pdf\n', outBase);
end

function make_fig(S, outBase)
W = 12; H = 8.5; left = 1.45; right = 0.1; top = 0.95; bottom = 0.95; gap = 1.9; legGap = 0.25;
lw = 1; lwRec = 2.2;
cTh = [0.8500 0.3250 0.0980]; cPhi = [0 0.4470 0.7410];
cThR = [0.4660 0.6740 0.1880]; cPhiR = [0.9290 0.6940 0.1250];
fig = figure('Units','centimeters','Position',[2 2 W H], 'Color','w', 'Visible', 'off');
set(fig, 'PaperUnits','centimeters', 'PaperSize',[W H], 'PaperPosition',[0 0 W H]);
boxH = (H - top - bottom - gap) / 2;
y1 = bottom + boxH + gap; y2 = bottom;
posRow = @(y, w) [left/W, y/H, w/W, boxH/H];
w0 = W - left - right - 4;
[b1, t1] = dual_x_row(posRow(y1, w0));
[b2, t2] = dual_x_row(posRow(y2, w0));
hAzT = plot(b1, S.phi, S.azT, 'Linewidth', lw, 'color', cPhi);
hAzR = plot(b1, S.phi, S.azR, '--', 'Linewidth', lwRec, 'color', cPhiR);
hElT = plot(t1, S.theta, S.elT, 'Linewidth', lw, 'color', cTh);
hElR = plot(t1, S.theta, S.elR, '--', 'Linewidth', lwRec, 'color', cThR);
ylabel(b1, S.ylab1, 'Interpreter','latex')
hAzE = plot(b2, S.phi, S.azE, 'Linewidth', lw, 'color', cPhi);
hElE = plot(t2, S.theta, S.elE, 'Linewidth', lw, 'color', cTh);
ylabel(b2, S.ylab2, 'Interpreter','latex')
for ax = [b1 b2], xlim(ax, S.phiLim); xlabel(ax, '$\phi$ ($^\circ$)', 'Interpreter','latex'); end
for ax = [t1 t2], xlim(ax, S.thetaLim); xlabel(ax, '$\theta$ ($^\circ$)', 'Interpreter','latex'); end
sync_ylim(b1, t1); sync_ylim(b2, t2);
lg1 = legend(b1, [hElT hElR hAzT hAzR], {'Truth $\theta$','Recovered $\theta$','Truth $\phi$','Recovered $\phi$'}, 'Interpreter','latex');
lg2 = legend(b2, [hElE hAzE], {'Error $\theta$','Error $\phi$'}, 'Interpreter','latex');
drawnow
lg1.Units = 'centimeters'; lg2.Units = 'centimeters';
legW = max(lg1.Position(3), lg2.Position(3));
boxW = W - left - right - legW - legGap;
for ax = [b1 t1], ax.Position = posRow(y1, boxW); end
for ax = [b2 t2], ax.Position = posRow(y2, boxW); end
lg1.Position(1:2) = [left + boxW + legGap, y1 + boxH - lg1.Position(4)];
lg2.Position(1:2) = [left + boxW + legGap, y2 + boxH - lg2.Position(4)];
saveas(fig, outBase, 'pdf'); close(fig);
fprintf('wrote %s.pdf\n', outBase);
end

function [axBot, axTop] = dual_x_row(pos)
axBot = axes('Position', pos);
axTop = axes('Position', pos, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none');
axTop.YTick = []; hold(axBot, 'on'); hold(axTop, 'on');
end

function sync_ylim(a, b)
yl = [min(a.YLim(1), b.YLim(1)), max(a.YLim(2), b.YLim(2))]; a.YLim = yl; b.YLim = yl;
end

function more_yticks(ax)
% Round tick step giving 4 to 6 labels.  Keeps MATLAB's limits when that works; otherwise extends
% a limit to the next multiple of the step, but only by less than a third of a step.
drawnow; yl = ax.YLim;
for s = [1 2 5 10 20 50 100 0.5 4]          % roundest steps first; 0.5 and 4 only if needed
    t = ceil(yl(1)/s - 1e-9)*s : s : yl(2) + 1e-9;
    if numel(t) >= 4 && numel(t) <= 6, ax.YTick = t; return; end
    lo = floor(yl(1)/s)*s; hi = ceil(yl(2)/s)*s;
    n = round((hi - lo)/s) + 1;
    if n >= 4 && n <= 6 && yl(1) - lo < s/3 && hi - yl(2) < s/3
        ax.YLim = [lo hi]; ax.YTick = lo:s:hi; return;
    end
end
end
