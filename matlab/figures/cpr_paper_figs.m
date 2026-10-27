function cpr_paper_figs(part)
%CPR_PAPER_FIGS  Section VI figures (cpr_fig.m: Times New Roman, LineWidth 1.5, no markers). Writes
%<results>/figs/paper/<name>.pdf/.png. Multi-panel figures are separate PDFs placed side by side in LaTeX.
%   cpr_paper_figs('IV')  pr2_ratio, pr2_P, pr2_K, pr2_rn, pr2_nd, pr2_snr, pr2_calerr (paper files pr3_*.pdf)
%   cpr_paper_figs('v9')  pr3_basiscmp (Fig. 6), pr3_calerr (Fig. 5(c)); each drawn only if its study has data
if nargin < 1, part = 'IV'; end
pDir = fullfile(cpr_paths().results, 'figs', 'paper');
if ~exist(pDir, 'dir'), mkdir(pDir); end
F = cpr_fig();
P = @(name) fullfile(pDir, name);
switch part
    case 'IV'
        % Section V panels: shorter (12 x 4.5 cm) for compact 2 x 2 figures; 11 pt text
        Fs = F; Fs.new = @() F.new(12, 4.5); Fs.save = @(f, file) F.save(f, file, 12, 4.5);
        % Fig. 6 (one row of three at 0.32\textwidth; the P panel is in the supplement):
        % drawn at 7.8 x 5 cm like the Fig. 5 panels, so the 11 pt text prints at full size
        F3 = F; F3.new = @() F.new(7.8, 5); F3.save = @(f, file) F.save(f, file, 7.8, 5);
        T = tryload('pr_ratio'); if ~isempty(T), archs(F3, T, 'ratio', 'metric', P('pr2_ratio'), 'sampling ratio {\itm/n}', 'northeast', true, []); end
        if ~isempty(T), archs(F3, T(ismember(T.arch, {'IS', 'RS'}), :), 'ratio', 'metric', P('pr2_ratio_ISRS'), 'sampling ratio {\itm/n}', 'northeast', true, []); end
        T = tryload('pr_P');     if ~isempty(T), single(F3, T, 'P', 'metric', P('pr2_P'), 'active ports {\itP}', [1 2 4 8 16], []); end
        T = tryload('pr_K');     if ~isempty(T), tones(F3, T, 'metric', P('pr2_K')); end
        T = tryload('pr_rn');    if ~isempty(T), rn(Fs, T, 'metric', P('pr2_rn'), []); end
        T = tryload('pr_nd');    if ~isempty(T), nd(F3, T, P('pr2_nd'), [-46 -20]); end
        T = tryload('pr_snr');   T2 = tryload('pr_snrseeds'); if ~isempty(T) && ~isempty(T2), T = [T; T2]; end   % IC/RC seeds 2-3 (mean of 3)
        if ~isempty(T), archs(F3, T, 'snrdB', 'metric', P('pr2_snr'), 'SNR (dB)', 'northeast', false, []); end
        T = tryload('pr_calerr'); if ~isempty(T), calerr(F3, T, P('pr2_calerr')); end
    case 'v9'
        % round 7 (2026-10-20): Fig. 7 (basis comparison, pr_basiscmp32) and Fig. 6(c) with self-calibration
        % (pr_calerr30 = the 30-update control without it, pr_calerrsc, pr_calerrics, pr_calerrrc), 7.8 x 5 cm like
        % the other Fig. 5 panels. Reads the current results folder (CPR_RESULTS).
        F3 = F; F3.new = @() F.new(7.8, 5); F3.save = @(f, file) F.save(f, file, 7.8, 5);
        T = tryload('pr_basiscmp32');
        if ~isempty(T), basiscmp(F, T, P('pr3_basiscmp')); end   % full-column 12 x 6 cm, like Fig. 10 (2026-10-21)
        Tsc = [tryload('pr_calerr30'); tryload('pr_calerrsc'); tryload('pr_calerrics'); tryload('pr_calerrrc')];
        if ~isempty(Tsc), calerrsc(F3, Tsc, P('pr3_calerr')); end
end
fprintf('cpr_paper_figs(%s): figures in %s\n', part, pDir);
end

% ---------------------------------------------------------------------------------------------
function m = metric(T, which)
% complex: raw; phaseless: per-wavelength alignment for K = 1, one global phase for K > 1
if strcmp(which, 'raw'), m = T.raw; return, end
m = T.global; k1 = T.K == 1; m(k1) = T.perlambda(k1);
end

function archs(F, T, xf, which, file, xl, loc, logx, yl)
f = F.new(); ax = axes(f); hold(ax, 'on'); h = []; nm = {};
for a = {'IS', 'RS', 'IC', 'RC'}
    s = T(strcmp(T.arch, a{1}), :); if isempty(s), continue; end
    s = sortrows(s, xf); y = metric(s, which);
    [xu, ~, g] = unique(s.(xf)); y0 = y; y = accumarray(g, y0, [], @mean);   % several seeds at one x: mean in dB
    h(end+1) = plot(ax, xu, y, F.ls.(a{1}), 'Color', F.col.(a{1}), 'LineWidth', F.lw); nm{end+1} = a{1}; %#ok<AGROW>
    if any(accumarray(g, 1) > 1)                                   % range over the trials at each x
        lo = accumarray(g, y0, [], @min); hi = accumarray(g, y0, [], @max);
        for i = 1:numel(xu), plot(ax, [xu(i) xu(i)], [lo(i) hi(i)], '-', 'Color', F.col.(a{1}), 'LineWidth', 0.75); end
    end
end
if logx
    set(ax, 'XScale', 'log', 'XTick', [0.03 0.05 0.1 0.2 0.4 0.8]);
    if ~F.paperAxes(), set(ax, 'XTick', unique(T.(xf))); end
end
xlim(ax, [min(T.(xf)) max(T.(xf))]); if ~isempty(yl) && F.paperAxes(), ylim(ax, yl); end
xlabel(ax, xl); ylabel(ax, 'NMSE (dB)'); if ~F.paperAxes(), loc = 'best'; end
legend(ax, h, nm, 'Location', loc); F.style(ax);
F.save(f, file);
end

function single(F, T, xf, which, file, xl, xt, yl)
T = sortrows(T, xf); f = F.new(); ax = axes(f);
plot(ax, T.(xf), metric(T, which), '-', 'Color', F.col.IS, 'LineWidth', F.lw);
if ~F.paperAxes(), xt = unique(T.(xf)).'; end
set(ax, 'XScale', 'log', 'XTick', xt); xlim(ax, [min(xt) max(xt)]); if ~isempty(yl) && F.paperAxes(), ylim(ax, yl); end
xlabel(ax, xl); ylabel(ax, 'NMSE (dB)'); F.style(ax); F.save(f, file);
end

function tones(F, T, which, file)
% NMSE only: the former right axis showed readings R = m/K, which follows
% directly from m/n and K and is stated in the text.
T = sortrows(T, 'K'); f = F.new(); ax = axes(f);
plot(ax, T.K, metric(T, which), '-', 'Color', F.col.IS, 'LineWidth', F.lw);
ylabel(ax, 'NMSE (dB)');
xt = [1 2 4 8 16]; if ~F.paperAxes(), xt = unique(T.K).'; end
set(ax, 'XScale', 'log', 'XTick', xt); xlim(ax, [min(xt) max(xt)]); xlabel(ax, 'tones per reading {\itK}');
F.style(ax); F.save(f, file);
end

function rn(F, T, which, file, yl)
f = F.new(); ax = axes(f); hold(ax, 'on'); h = []; nm = {}; xall = [];
for a = {'RS', 'RC'}
    s = T(strcmp(T.arch, a{1}), :); if isempty(s), continue; end
    xv = s.rnRS; if strcmp(a{1}, 'RC'), xv = s.rnRC; end
    [xv, o] = sort(xv); y = metric(s, which); y = y(o); xall = [xall; xv(:)]; %#ok<AGROW>
    h(end+1) = plot(ax, xv, y, F.ls.(a{1}), 'Color', F.col.(a{1}), 'LineWidth', F.lw); nm{end+1} = a{1}; %#ok<AGROW>
end
xt = [1/32 1/16 1/8 1/4 1/2 1]; if ~F.paperAxes(), xt = unique(xall).'; end
set(ax, 'XScale', 'log', 'XTick', xt, 'XTickLabel', arrayfun(@(v) fracLabel(v), xt, 'UniformOutput', false));
xlim(ax, [min(xt) max(xt)]); if ~isempty(yl) && F.paperAxes(), ylim(ax, yl); end
xlabel(ax, 'installed-probe fraction {\itr_n}'); ylabel(ax, 'NMSE (dB)'); legend(ax, h, nm, 'Location', ifelse(F.paperAxes(), 'northeast', 'best'));
F.style(ax); F.save(f, file);
end

function nd(F, T, file, yl)
which = 'metric'; if strcmp(T.meas{1}, 'complex'), which = 'raw'; end
j = startsWith(T.label, 'joint'); s1 = startsWith(T.label, 'single');
Tj = sortrows(T(j, :), 'nd'); Ts = T(s1, :); fi = cellfun(@(l) sscanf(l, 'single_f%d'), Ts.label);
ys = metric(Ts, which);
% "each alone": the pooled NMSE of the independently recovered wavelengths, i.e. their errors summed
% and divided by the summed signal energy, the same aggregation as the joint NMSE
c0 = cpr_cfg();
if strcmp(c0.angBasis, 'vsh')                      % the recovered (exact, possibly dual-polarized) pattern
    Xe = pattern_crlh_vsh(Ts.nw(1), max(fi), 1, 16, c0.npol);
    Ek = squeeze(sum(sum(abs(Xe).^2, 1), 2)); Ek = sum(reshape(Ek, max(fi), c0.npol), 2);
else
    Xe = pattern_crlh(Ts.nw(1), max(fi), 1); Ek = squeeze(sum(sum(abs(Xe).^2, 1), 2)); Ek = Ek(:);
end
yi = arrayfun(@(n) 10*log10(sum(10.^(ys(fi <= n)/10) .* Ek(fi(fi <= n))) / sum(Ek(fi(fi <= n)))), Tj.nd);
f = F.new(); ax = axes(f); hold(ax, 'on');
plot(ax, Tj.nd, metric(Tj, which), '-', 'Color', F.col.IS, 'LineWidth', F.lw);
plot(ax, Tj.nd, yi, '--', 'Color', F.col.RS, 'LineWidth', F.lw);
xt = Tj.nd.'; set(ax, 'XScale', 'log', 'XTick', xt); xlim(ax, [min(xt) max(xt)]); if ~isempty(yl) && F.paperAxes(), ylim(ax, yl); end
xlabel(ax, 'wavelengths {\itn}_\lambda'); ylabel(ax, 'NMSE (dB)');
loc = 'northeast'; if strcmp(which, 'metric'), loc = 'north'; end   % phaseless: both curves start top left and "each alone" ends top right
legend(ax, {'joint', 'each alone'}, 'Location', loc); F.style(ax); F.save(f, file);
end

function calerr(F, T, file)
% NMSE vs rms cable phase error sigma_phi for IS and IC, with the 20 log10(sigma_phi) reference
f = F.new(); ax = axes(f); hold(ax, 'on'); h = []; nm = {};
for a = {'IS', 'IC'}
    s = sortrows(T(strcmp(T.arch, a{1}), :), 'calErrDeg');
    [xu, ~, g] = unique(s.calErrDeg); y = accumarray(g, s.global, [], @mean);
    h(end+1) = plot(ax, xu, y, F.ls.(a{1}), 'Color', F.col.(a{1}), 'LineWidth', F.lw); nm{end+1} = a{1}; %#ok<AGROW>
end
xr = linspace(0.5, max(T.calErrDeg), 50);
h(end+1) = plot(ax, xr, 20*log10(xr*pi/180), ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
nm{end+1} = '20log_{10}\sigma_\phi';
xlim(ax, [0 max(T.calErrDeg)]);
xlabel(ax, 'rms cable phase error \sigma_\phi (deg)'); ylabel(ax, 'NMSE (dB)');
legend(ax, h, nm, 'Location', 'southeast'); F.style(ax); F.save(f, file);
end

function basiscmp(F, T, file)
% NMSE vs sampling ratio: colour = architecture, solid = spherical-wave model, dashed = unconstrained 2D Fourier model
f = F.new(); ax = axes(f); hold(ax, 'on'); h = []; nm = {};
% IS and IC coincide with the spherical-wave model: IC is drawn dotted and thicker on top of the solid IS line, so both
% stay visible where they overlap (user, 2026-10-21); the Fourier curves are dashed (IS) and dash-dot (IC)
cur = {'IS_vsh', 'IS', '-', 'IS, spherical-wave', F.lw; 'IS_dft', 'IS', '--', 'IS, Fourier', F.lw; ...
       'IC_vsh', 'IC', ':', 'IC, spherical-wave', 2.2; 'IC_dft', 'IC', '-.', 'IC, Fourier', F.lw};
for i = 1:size(cur, 1)
    s = T(startsWith(T.label, cur{i, 1}), :); if isempty(s), continue; end
    s = sortrows(s, 'ratio');
    [xu, ~, g] = unique(s.ratio); y = accumarray(g, metric(s, 'metric'), [], @mean);
    h(end+1) = plot(ax, xu, y, cur{i, 3}, 'Color', F.col.(cur{i, 2}), 'LineWidth', cur{i, 5}); nm{end+1} = cur{i, 4}; %#ok<AGROW>
end
set(ax, 'XScale', 'log', 'XTick', [0.05 0.1 0.2 0.4 0.8]);
if ~F.paperAxes(), set(ax, 'XTick', unique(T.ratio)); end
xlim(ax, [min(T.ratio) max(T.ratio)]);
xlabel(ax, 'sampling ratio {\itm/n}'); ylabel(ax, 'NMSE (dB)');
F.style(ax);
yl = ylim(ax); top = 10*ceil(yl(2)/10); bot = 20*floor(yl(1)/20);   % headroom above the curves (top 10 dB here)
ylim(ax, [bot - 40, top]); ax.YTick = bot:20:0;                     % empty strip at the bottom for the legend
legend(ax, h, nm, 'Location', 'south', 'NumColumns', 2); F.save(f, file);
end

function calerrsc(F, T, file)
% NMSE vs rms cable phase error at the band centre: with self-calibration (solid) for IS, IC, two-pass IC (label
% prefix ICS) and RC; the equal-effort control without it (dashed) for IS and IC; and the expected error of a
% conventional measurement through the same cables (dotted), E|exp(i e_k) - 1|^2 averaged over the 16 frequencies
% with e_k ~ N(0, (sigma f_k / f_c)^2) (length errors: the phase grows with frequency).
f = F.new(); ax = axes(f); hold(ax, 'on'); h = []; nm = {};
grp = regexp(T.label, '^[A-Z]+', 'match', 'once');           % IS, IC, ICS, RC
sc = endsWith(T.label, '_sc');
% Legend: the four colours only (solid = self-calibration, dashed = without, dotted = conventional; stated in the
% caption), two columns in a strip below the curves, so it covers no data (user, 2026-10-21).
cur = {'IS', true, 'IS', '-', 'IS'; 'IC', true, 'IC', '-', 'IC'; ...
       'ICS', true, 'RS', '-', 'two-pass IC'; 'RC', true, 'RC', '-', 'RC'; ...
       'IS', false, 'IS', '--', ''; 'IC', false, 'IC', '--', ''};
for i = 1:size(cur, 1)
    s = T(strcmp(grp, cur{i, 1}) & sc == cur{i, 2} & T.calErrDeg > 0, :); if isempty(s), continue; end
    [xu, ~, g] = unique(s.calErrDeg); y = accumarray(g, s.global, [], @mean);
    hi = plot(ax, xu, y, cur{i, 4}, 'Color', F.col.(cur{i, 3}), 'LineWidth', F.lw);
    if ~isempty(cur{i, 5}), h(end+1) = hi; nm{end+1} = cur{i, 5}; end %#ok<AGROW>
end
fk = linspace(1.90, 2.65, 16); fc = mean(fk);
xr = linspace(0.5, max(T.calErrDeg), 60);
yc = arrayfun(@(sd) 10*log10(mean(2*(1 - exp(-((sd*pi/180)*fk/fc).^2/2)))), xr);
plot(ax, xr, yc, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
xlim(ax, [0 max(T.calErrDeg)]);
xlabel(ax, 'rms cable phase error \sigma_\phi (deg)'); ylabel(ax, 'NMSE (dB)');
F.style(ax);
yl = ylim(ax); ylim(ax, [yl(1) - 15, yl(2)]); ax.YTick = ceil(yl(1)/10)*10 : 10 : floor(yl(2)/10)*10;
legend(ax, h, nm, 'Location', 'south', 'NumColumns', 2); F.save(f, file);
end

function T = tryload(study)
% studies that have not produced results yet are skipped
try, T = cpr_load(study); catch, T = []; end
end

function s = fracLabel(v)
% 1 or 1/N for the installed-probe fractions
if abs(v - 1) < 1e-9, s = '1'; else, s = sprintf('1/%g', round(1/v)); end
end

function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end
