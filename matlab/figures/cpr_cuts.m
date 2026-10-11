function cpr_cuts(k, study, cutdir)
%CPR_CUTS  Truth vs. recovered pattern cut through the pattern peak, with the pointwise error, for
%the IS examples at two NMSE levels (paper Fig. 4 and supplementary Figs. S2 and S6).
%   cpr_cuts(1, 'pr')          phaseless examples (pr_examples), elevation cut at the first frequency
%   cpr_cuts(1, 'pr', 'phi')   the same examples along the azimuth cut through the peak
%   cpr_cuts(1, 'cs')          complex-measurement examples (cs_examples)
%cutdir 'theta' (default) is the elevation cut through the peak. 'phi' is the azimuth cut at the
%peak's polar angle: the scan grid spans phi 0-180 deg, so the scan rows theta and 360 - theta
%are joined to give phi 0-360 deg (output names get '_phi'). The two examples are the runs whose
%whole-pattern NMSE is nearest -20 and -34 dB (complex) or -24 and -36 dB (phaseless). Layout: 12 cm
%wide, Times 11 pt, four stacked panels (magnitude, its error, phase, its error), legend at the right.
%Writes <results>/figs/paper/cs2_cut*.pdf or pr2_cut*.pdf.
if nargin < 3, cutdir = 'theta'; end
outDir = cpr_paths().results;
pDir = fullfile(outDir, 'figs', 'paper');
set(0,'defaultAxesFontName','Times'); set(0,'defaultTextFontName','Times');
set(0,'defaultTextFontSize',11);      set(0,'defaultAxesFontSize',11);
set(0,'defaultLegendFontName','Times');
% study 'cs' (default): Figs. 6-7, complex IS examples at -20/-34 dB (raw NMSE; complex data have no
% phase ambiguity). 'pr': Section IV.A phaseless IS examples at -24/-36 dB (K = 4, so one global
% phase is removed over all wavelengths). Outputs cs2_cut* or pr2_cut*.
if nargin < 1 || isempty(k), k = 1; end   % 1.90 GHz: most detailed theta cut (cpr_pick_cut.m)
if nargin < 2, study = 'cs'; end
switch study
    case 'cs', targets = [-20 -34]; scoreF = 'raw';    pre = 'cs2';
    case 'pr', targets = [-24 -36]; scoreF = 'global'; pre = 'pr2';
end
E = {};
for d = dir(fullfile(outDir, [study '_examples_w*.mat'])).'
    S0 = load(fullfile(d.folder, d.name));
    for e = S0.res(:).', if ~isempty(e.Xh), E{end+1} = e; end, end %#ok<AGROW>
end
nm = cellfun(@(e) e.(scoreF), E);
[~, i1] = min(abs(nm - targets(1))); [~, i2] = min(abs(nm - targets(2))); E = E([i1 i2]);
c = E{1}.cfg; [X, ~] = pattern_crlh(c.nw, c.nd, c.f0idx);
for i = 1:2                                % remove one global phase over all wavelengths
    Xh = E{i}.Xh; E{i}.Xh = Xh * exp(1i*angle(Xh(:)'*X(:))); E{i}.raw = E{i}.(scoreF);
end
US = X(:, :, k); nw = size(US, 1);
[~, ip] = max(abs(US(:))); [~, slicep1] = ind2sub(size(US), ip);        % elevation cut through the peak
theta_deg = (0:nw-1)*360/nw;
peakMag = max(abs(US(:)));
errdB   = @(a, b) max(20*log10(abs(a - b) / peakMag), -100);
wrapdeg = @(a) angle(exp(1i*a)) * 180/pi;
mn = char(8722);
S = struct(); S.theta = theta_deg; S.thetaLim = [0 360-360/nw];
S.nm = {sprintf('%s%.0f dB', mn, -E{1}.raw), sprintf('%s%.0f dB', mn, -E{2}.raw)};
elT = US(:, slicep1);
for i = 1:2, d = E{i}.Xh(:,:,k); d = d*exp(1i*angle(d(:)'*US(:))); fprintf('example %d (NMSE %.1f dB overall): %.1f dB at %.2f GHz, theta cut at phi = %.2f deg\n', i, E{i}.raw, 20*log10(norm(d(:)-US(:))/norm(US(:))), 1.9 + 0.05*(k-1), (slicep1-1)*180/c.nh); end
elR = cellfun(@(e) double(e.Xh(:, slicep1, k)), E, 'UniformOutput', false);
% --- magnitude (Fig. 6) ---
S.elT = 20*log10(abs(elT)/peakMag);
S.elR = cellfun(@(v) 20*log10(abs(v)/peakMag), elR, 'UniformOutput', false);
S.elE = cellfun(@(v) errdB(v, elT), elR, 'UniformOutput', false);
S.ylab1 = '$\left| E \right|$ (dB)';
S.ylab2 = '$\left| \mathrm{Error} \right|$ (dB)';
make_fig(S, fullfile(pDir, [pre '_cutmag']));
% --- phase (Fig. 7) ---
S.elT = angle(elT)*180/pi;
S.elR = cellfun(@(v) angle(v)*180/pi, elR, 'UniformOutput', false);
S.elE = cellfun(@(v) wrapdeg(angle(v) - angle(elT)), elR, 'UniformOutput', false);
S.ylab1 = '$\angle E$ ($^\circ$)';
S.ylab2 = 'Phase error ($^\circ$)';
make_fig(S, fullfile(pDir, [pre '_cutphase']));
% --- combined magnitude + phase, four stacked panels (one figure per pair) ---
C = S; C.xlab = '$\theta$ ($^\circ$)'; outName = [pre '_cut'];
if strcmp(cutdir, 'phi')                   % azimuth cut through the peak (supplement)
    [ith, ~] = ind2sub(size(US), ip); nh = size(US, 2);
    ith2 = mod(nw - (ith - 1), nw) + 1;    % scan theta' = 360 - theta: same polar angle, phi + 180
    if (ith - 1)*360/nw > 180, [ith, ith2] = deal(ith2, ith); end   % row with scan theta <= 180 gives phi 0-180
    cut = @(M) [M(ith, :), M(ith2, :)].';
    elT = cut(US); elR = cellfun(@(e) cut(double(e.Xh(:, :, k))), E, 'UniformOutput', false);
    C.theta = (0:2*nh-1)*180/nh; C.thetaLim = [0 360-180/nh]; C.xlab = '$\phi$ ($^\circ$)';
    outName = [pre '_cut_phi'];
    fprintf('azimuth cut at scan theta = %.2f deg (rows %d and %d)\n', (ith-1)*360/nw, ith, ith2);
end
C.elT = elT;
C.mT = 20*log10(abs(elT)/peakMag);
C.mR = cellfun(@(v) 20*log10(abs(v)/peakMag), elR, 'UniformOutput', false);
C.mE = cellfun(@(v) errdB(v, elT), elR, 'UniformOutput', false);
C.pT = angle(elT)*180/pi;
C.pR = cellfun(@(v) angle(v)*180/pi, elR, 'UniformOutput', false);
C.pE = cellfun(@(v) wrapdeg(angle(v) - angle(elT)), elR, 'UniformOutput', false);
make_fig4(C, fullfile(pDir, outName));
end

function make_fig4(S, outBase)
    W = 12; H = 10.5; left = 1.6; right = 0.1; top = 0.2; bottom = 1.0; gap = 0.3; legGap = 0.25;
    lw = 1; lwRec = 1.6;
    cTruth = [0 0.4470 0.7410]; cRec = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250]; ls = {'--', '-.'};
    fig = figure('Visible','off','Units','centimeters','Position',[2 2 W H], 'Color','w');
    set(fig, 'PaperUnits','centimeters', 'PaperSize',[W H], 'PaperPosition',[0 0 W H]);
    boxH = (H - top - bottom - 3*gap) / 4;
    y = @(k) bottom + (4 - k)*(boxH + gap);
    w0 = W - left - right - 3;
    ax = gobjects(1, 4);
    for k = 1:4, ax(k) = axes('Units','centimeters','Position',[left, y(k), w0, boxH]); hold(ax(k), 'on'); end
    hT = plot(ax(1), S.theta, S.mT, 'Linewidth', lw, 'color', cTruth);
    hP = plot(ax(3), S.theta, S.pT, 'Linewidth', lw, 'color', cTruth);
    hR = gobjects(1, 2);
    for i = 1:2
        hR(i) = plot(ax(1), S.theta, S.mR{i}, ls{i}, 'Linewidth', lwRec, 'color', cRec(i, :));
        plot(ax(2), S.theta, S.mE{i}, ls{i}, 'Linewidth', 1.2, 'color', cRec(i, :));
        plot(ax(3), S.theta, S.pR{i}, ls{i}, 'Linewidth', lwRec, 'color', cRec(i, :));
        plot(ax(4), S.theta, S.pE{i}, ls{i}, 'Linewidth', 1.2, 'color', cRec(i, :));
    end
    uistack(hT, 'top'); uistack(hP, 'top');           % truth drawn last in both panels
    ylabel(ax(1), '$\left| E \right|$ (dB)', 'Interpreter', 'latex');
    ylabel(ax(2), '$\left| \mathrm{Error} \right|$ (dB)', 'Interpreter', 'latex');
    ylabel(ax(3), '$\angle E$ ($^\circ$)', 'Interpreter', 'latex');
    ylabel(ax(4), 'Error ($^\circ$)', 'Interpreter', 'latex');   % short: fits the shorter panels
    ylim(ax(3), [-180 180]); ax(3).YTick = -180:90:180;
    for k = 1:4, xlim(ax(k), S.thetaLim); box(ax(k), 'on'); end
    for k = 1:3, ax(k).XTickLabel = []; end
    for k = [1 2 4], more_yticks(ax(k)); end         % 4-6 labels per panel
    xlabel(ax(4), S.xlab, 'Interpreter', 'latex');
    lg = legend(ax(1), [hT hR], [{'Truth'}, S.nm]); drawnow
    lg.Units = 'centimeters'; legW = 1.08*lg.Position(3);
    boxW = W - left - right - legW - legGap;
    for k = 1:4, ax(k).Position = [left, y(k), boxW, boxH]; end
    lg.Position(3) = legW; lg.Position(1:2) = [left + boxW + legGap, y(1) + boxH - lg.Position(4)];
    saveas(fig, outBase, 'pdf');
    fprintf('wrote %s.pdf\n', outBase);
    close(fig);
end

% =========================================================================
function make_fig(S, outBase)
    W = 12;  H = 7;                   % page size, cm (width matches paper figs)
    left = 1.45; right = 0.1;         % cm
    top = 0.25; bottom = 0.95;        % room for the theta label below each row
    gap = 1.15;                       % between rows: row 1 theta label
    legGap = 0.25;                    % cm between plot box and legend

    lw = 1;  lwRec = 1.6;
    cTruth = [0 0.4470 0.7410];       % truth: MATLAB color 1 (blue)
    cRec   = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];   % -20 dB red, -34 dB yellow (colors 2, 3),
                                                             % same for their errors
    ls = {'--', '-.'};                % -20 dB dashed, -34 dB dash-dot

    fig = figure('Visible','off','Units','centimeters','Position',[2 2 W H], 'Color','w');
    set(fig, 'PaperUnits','centimeters', 'PaperSize',[W H], 'PaperPosition',[0 0 W H]);

    boxH = (H - top - bottom - gap) / 2;
    y1 = bottom + boxH + gap;   y2 = bottom;
    posRow = @(y, w) [left/W, y/H, w/W, boxH/H];

    w0 = W - left - right - 3;
    a1 = axes('Position', posRow(y1, w0)); hold(a1, 'on');
    a2 = axes('Position', posRow(y2, w0)); hold(a2, 'on');

    hT = plot(a1, S.theta, S.elT, 'Linewidth', lw, 'color', cTruth);
    hR = gobjects(1, 2); hE = hR;
    for i = 1:2
        hR(i) = plot(a1, S.theta, S.elR{i}, ls{i}, 'Linewidth', lwRec, 'color', cRec(i, :));
        hE(i) = plot(a2, S.theta, S.elE{i}, ls{i}, 'Linewidth', 1.2, 'color', cRec(i, :));
    end
    uistack(hT, 'top');
    ylabel(a1, S.ylab1, 'Interpreter', 'latex'); ylabel(a2, S.ylab2, 'Interpreter', 'latex');   % axis labels: LaTeX as Figs. 18-19; legends stay plain TeX
    for ax = [a1 a2]
        xlim(ax, S.thetaLim); xlabel(ax, '$\theta$ ($^\circ$)', 'Interpreter', 'latex'); box(ax, 'on');
    end

    lg1 = legend(a1, [hT hR], [{'Truth'}, S.nm]);
    lg2 = legend(a2, hE, {[S.nm{1} ' error'], [S.nm{2} ' error']});
    drawnow

    lg1.Units = 'centimeters';  lg2.Units = 'centimeters';
    legW = 1.08 * max(lg1.Position(3), lg2.Position(3));     % small margin so the lines fit
    boxW = W - left - right - legW - legGap;
    a1.Position = posRow(y1, boxW); a2.Position = posRow(y2, boxW);
    lg1.Position(3) = legW; lg2.Position(3) = legW;
    lg1.Position(1:2) = [left + boxW + legGap, y1 + boxH - lg1.Position(4)];
    lg2.Position(1:2) = [left + boxW + legGap, y2 + boxH - lg2.Position(4)];

    saveas(fig, outBase, 'pdf')
    print(fig, [outBase '.png'], '-dpng', '-r300')
    fprintf('wrote %s.pdf / .png\n', outBase);
    close(fig);
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
