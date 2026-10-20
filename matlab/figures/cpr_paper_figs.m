function cpr_paper_figs(part)
%CPR_PAPER_FIGS  Paper figures (cpr_fig.m: Times New Roman, LineWidth 1.5, no markers). Writes
%<results>/figs/paper/<name>.pdf/.png. Multi-panel figures are separate PDFs placed side by side in
%LaTeX.
%   cpr_paper_figs('supplement')  supplementary complex measurements:  cs2_A (24 x 6 cm), cs2_sparsity, cs2_ratio,
%                         cs2_P, cs2_K, cs2_rn, cs2_nd, cs2_snr, cs2_size, cs2_exmag, cs2_exphase
%   cpr_paper_figs('IV')   Section V (phaseless): pr2_ratio, pr2_P, pr2_K, pr2_rn, pr2_nd,
%                         pr2_legs, pr2_snr, pr2_size (each drawn only if its study has data)
if nargin < 1, part = 'supplement'; end
pDir = fullfile(cpr_paths().results, 'figs', 'paper');
if ~exist(pDir, 'dir'), mkdir(pDir); end
F = cpr_fig();
P = @(name) fullfile(pDir, name);
switch part
    case 'supplement'
        matrices(F, P('cs2_A'));
        sparsity(F, P('cs2_sparsity'));
        archs(F, cpr_load('cs_ratio'), 'ratio', 'raw', P('cs2_ratio'), 'sampling ratio {\itm/n}', 'northeast', true, []);
        % Fig. 6 panels: 3 per row at 0.32\textwidth (5.8 cm), drawn at
        % 7.8 x 5 cm so the 11 pt text prints at the same size as the 12 cm figures at 0.49\textwidth.
        Fs = F; Fs.new = @() F.new(7.8, 5); Fs.save = @(f, file) F.save(f, file, 7.8, 5);
        single(Fs, cpr_load('cs_P'), 'P', 'raw', P('cs2_P'), 'active probes {\itP}', [1 2 4 8 16 32], [-50 -10]);
        tones(Fs, cpr_load('cs_K'), 'raw', P('cs2_K'));
        rn(Fs, cpr_load('cs_rn'), 'raw', P('cs2_rn'), [-50 0]);
        nd(Fs, cpr_load('cs_nd'), P('cs2_nd'), [-40 -10]);
        archs(Fs, cpr_load('cs_snr'), 'snrdB', 'raw', P('cs2_snr'), 'SNR (dB)', 'northeast', false, [-40 -10]);
        single(Fs, cpr_load('cs_size'), 'nw', 'raw', P('cs2_size'), '{\itn}_\theta  ({\itn}_\phi = {\itn}_\theta/2)', [16 32 64 128], [-50 0]);
        examples(F, 'cs_examples', P('cs2_exmag'), P('cs2_exphase'));
    case 'supplementsweeps'                              % only the condensed Fig. 6 panels (same as in 'supplement')
        Fs = F; Fs.new = @() F.new(7.8, 5); Fs.save = @(f, file) F.save(f, file, 7.8, 5);
        single(Fs, cpr_load('cs_P'), 'P', 'raw', P('cs2_P'), 'active probes {\itP}', [1 2 4 8 16 32], [-50 -10]);
        tones(Fs, cpr_load('cs_K'), 'raw', P('cs2_K'));
        rn(Fs, cpr_load('cs_rn'), 'raw', P('cs2_rn'), [-50 0]);
        nd(Fs, cpr_load('cs_nd'), P('cs2_nd'), [-40 -10]);
        archs(Fs, cpr_load('cs_snr'), 'snrdB', 'raw', P('cs2_snr'), 'SNR (dB)', 'northeast', false, [-40 -10]);
        single(Fs, cpr_load('cs_size'), 'nw', 'raw', P('cs2_size'), '{\itn}_\theta  ({\itn}_\phi = {\itn}_\theta/2)', [16 32 64 128], [-50 0]);
    case 'fig6'                                    % only the NMSE-example figure (Fig. 6)
        examples(F, 'cs_examples', P('cs2_exmag'), P('cs2_exphase'));
    case 'IVnd'                                    % only the phaseless n_lambda panel
        F3 = F; F3.new = @() F.new(7.8, 5); F3.save = @(f, file) F.save(f, file, 7.8, 5);
        nd(F3, cpr_load('pr_nd'), P('pr2_nd'), [-46 -20]);
    case 'IVP'                                     % only the phaseless P panel
        F3 = F; F3.new = @() F.new(7.8, 5); F3.save = @(f, file) F.save(f, file, 7.8, 5);
        single(F3, cpr_load('pr_P'), 'P', 'metric', P('pr2_P'), 'active ports {\itP}', [1 2 4 8 16], []);
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
        T = tryload('pr_legs');  if ~isempty(T), archs(F3, T, 'nlegs', 'metric', P('pr2_legs'), 'number of leg probes', 'northeast', false, []); end
        T = tryload('pr_snr');   T2 = tryload('pr_snrseeds'); if ~isempty(T) && ~isempty(T2), T = [T; T2]; end   % IC/RC seeds 2-3 (mean of 3)
        if ~isempty(T), archs(F3, T, 'snrdB', 'metric', P('pr2_snr'), 'SNR (dB)', 'northeast', false, []); end
        T = tryload('pr_calerr'); if ~isempty(T), calerr(F3, T, P('pr2_calerr')); end
        T = tryload('pr_size');  if ~isempty(T), single(F3, T, 'nw', 'metric', P('pr2_size'), '{\itn}_\theta  ({\itn}_\phi = {\itn}_\theta/2)', [16 32 64 128], []); end
    case 'v9'
        % round 7 (2026-10-20): Fig. 7 (basis comparison, pr_basiscmp32) and Fig. 6(c) with self-calibration
        % (pr_calerr30 = the 30-update control without it, pr_calerrsc, pr_calerrics, pr_calerrrc), 7.8 x 5 cm like
        % the other Fig. 6 panels. Reads the current results folder (CPR_RESULTS); pr_basiscmp32 may live in another
        % folder: cpr_paper_figs('v9') looks for it in the current folder first, then in a sibling results_v8.
        F3 = F; F3.new = @() F.new(7.8, 5); F3.save = @(f, file) F.save(f, file, 7.8, 5);
        T = tryload('pr_basiscmp32');
        if isempty(T)
            alt = fullfile(fileparts(cpr_paths().results), 'results_v8');
            if exist(alt, 'dir'), T = load_from(alt, 'pr_basiscmp32'); end
        end
        if ~isempty(T), basiscmp(F3, T, P('pr3_basiscmp')); end
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
cur = {'IS_vsh', 'IS', '-', 'IS, spherical-wave'; 'IS_dft', 'IS', '--', 'IS, Fourier'; ...
       'IC_vsh', 'IC', '-', 'IC, spherical-wave'; 'IC_dft', 'IC', '--', 'IC, Fourier'};
for i = 1:size(cur, 1)
    s = T(startsWith(T.label, cur{i, 1}), :); if isempty(s), continue; end
    s = sortrows(s, 'ratio');
    [xu, ~, g] = unique(s.ratio); y = accumarray(g, metric(s, 'metric'), [], @mean);
    h(end+1) = plot(ax, xu, y, cur{i, 3}, 'Color', F.col.(cur{i, 2}), 'LineWidth', F.lw); nm{end+1} = cur{i, 4}; %#ok<AGROW>
end
set(ax, 'XScale', 'log', 'XTick', [0.05 0.1 0.2 0.4 0.8]);
if ~F.paperAxes(), set(ax, 'XTick', unique(T.ratio)); end
xlim(ax, [min(T.ratio) max(T.ratio)]);
xlabel(ax, 'sampling ratio {\itm/n}'); ylabel(ax, 'NMSE (dB)');
legend(ax, h, nm, 'Location', ifelse(F.paperAxes(), 'southwest', 'best')); F.style(ax); F.save(f, file);
end

function calerrsc(F, T, file)
% NMSE vs rms cable phase error at the band centre: with self-calibration (solid) for IS, IC, two-pass IC (label
% prefix ICS) and RC; the equal-effort control without it (dashed) for IS and IC; and the expected error of a
% conventional measurement through the same cables (dotted), E|exp(i e_k) - 1|^2 averaged over the 16 frequencies
% with e_k ~ N(0, (sigma f_k / f_c)^2) (length errors: the phase grows with frequency).
f = F.new(); ax = axes(f); hold(ax, 'on'); h = []; nm = {};
grp = regexp(T.label, '^[A-Z]+', 'match', 'once');           % IS, IC, ICS, RC
sc = endsWith(T.label, '_sc');
cur = {'IS', true, 'IS', '-', 'IS, self-cal.'; 'IC', true, 'IC', '-', 'IC, self-cal.'; ...
       'ICS', true, 'RS', '-', 'two-pass IC, self-cal.'; 'RC', true, 'RC', '-', 'RC, self-cal.'; ...
       'IS', false, 'IS', '--', 'IS'; 'IC', false, 'IC', '--', 'IC'};
for i = 1:size(cur, 1)
    s = T(strcmp(grp, cur{i, 1}) & sc == cur{i, 2} & T.calErrDeg > 0, :); if isempty(s), continue; end
    [xu, ~, g] = unique(s.calErrDeg); y = accumarray(g, s.global, [], @mean);
    h(end+1) = plot(ax, xu, y, cur{i, 4}, 'Color', F.col.(cur{i, 3}), 'LineWidth', F.lw); nm{end+1} = cur{i, 5}; %#ok<AGROW>
end
fk = linspace(1.90, 2.65, 16); fc = mean(fk);
xr = linspace(0.5, max(T.calErrDeg), 60);
yc = arrayfun(@(sd) 10*log10(mean(2*(1 - exp(-((sd*pi/180)*fk/fc).^2/2)))), xr);
h(end+1) = plot(ax, xr, yc, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1); nm{end+1} = 'conventional';
xlim(ax, [0 max(T.calErrDeg)]);
xlabel(ax, 'rms cable phase error \sigma_\phi (deg)'); ylabel(ax, 'NMSE (dB)');
legend(ax, h, nm, 'Location', 'southeast', 'FontSize', 7); F.style(ax); F.save(f, file);
end

function T = load_from(dirName, study)
% cpr_load from another results folder
old = getenv('CPR_RESULTS'); setenv('CPR_RESULTS', dirName);
try, T = cpr_load(study); catch, T = []; end
setenv('CPR_RESULTS', old);
end

function matrices(F, file)
D = crlh_data();
fq = D.freq(1:2).';
f = F.new(24, 6); tl = tiledlayout(f, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
a4 = {'IS', 'RS', 'IC', 'RC'}; ttl = {'(a) IS', '(b) RS', '(c) IC', '(d) RC'};
for i = 1:4
    c = cpr_cfg('meas', 'complex', 'arch', a4{i}, 'nw', 8, 'nd', 2, 'K', 2, 'P', 2, 'ratio', 0.25, ...
        'nlegs', 0, 'rnRS', 1/2, 'rnRC', 1/2);
    S = make_readings(c, fq, 3); n2 = c.nw*c.nh;
    ax = nexttile(tl); imagesc(ax, ~full(S.Phi ~= 0)); colormap(ax, gray); axis(ax, 'normal'); hold(ax, 'on');
    plot(ax, [n2 n2] + 0.5, [0.5 S.m + 0.5], 'r-', 'LineWidth', 1);
    for q = 1:c.nh - 1
        plot(ax, [q q]*c.nw + 0.5, [0.5 S.m + 0.5], ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
        plot(ax, ([q q] + c.nh)*c.nw + 0.5, [0.5 S.m + 0.5], ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
    end
    for r = 2:2:S.m - 1, plot(ax, [0.5 2*n2 + 0.5], [r r] + 0.5, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.25); end
    set(ax, 'XTick', [1 n2 + 1], 'XTickLabel', {'\lambda_1', '\lambda_2'}, 'YTick', [], 'FontName', F.font);
    title(ax, ttl{i}, 'FontWeight', 'normal', 'FontName', F.font);
    if i == 1, ylabel(ax, 'measurements'); end
end
xlabel(tl, 'pattern samples', 'FontName', F.font);
F.save(f, file, 24, 6);
end

function sparsity(F, file)
c = cpr_cfg(); [X, fr] = pattern_crlh(c.nw, c.nd, c.f0idx);
[nw, nh, nd_] = size(X); n2 = nw*nh; Xm = reshape(X, n2, nd_);
S2 = reshape(fft2(reshape(Xm, nw, nh, [])), n2, []) / sqrt(n2);
rng(1); rp = exp(2i*pi*rand(1, nd_));
dct = cpr_basis('dct', fr, 0);
tau = -angle(sum(sum(conj(S2(:, 1:end-1)) .* S2(:, 2:end)))) / (2*pi*(fr(2) - fr(1)));
Mt = cpr_basis('dct', fr, tau);
cv = {S2, fft(S2, [], 2)/sqrt(nd_), fft(S2 .* rp, [], 2)/sqrt(nd_), (dct' * S2.').', (Mt' * S2.').'};
nm = {'each wavelength alone', '3D DFT', '3D DFT, random phase', '2D DFT \times DCT', '2D DFT \times DCT, delay comp.'};
f = F.new(); ax = axes(f); hold(ax, 'on');
for k = 1:5
    e = sort(abs(cv{k}(:)).^2, 'descend'); ce = cumsum(e)/sum(e);
    plot(ax, 1:numel(ce), 1 - ce, F.lss{k}, 'Color', F.cols(k, :), 'LineWidth', F.lw);
end
set(ax, 'XScale', 'log', 'YScale', 'log', 'YTick', 10.^(-5:0), 'XTick', 10.^(0:4)); ylim(ax, [1e-5 1]); xlim(ax, [1 1e4]);
xlabel(ax, 'largest coefficients kept'); ylabel(ax, 'energy not captured');
legend(ax, nm, 'Location', 'southwest', 'FontSize', 7); F.style(ax); F.save(f, file);
end

function examples(F, study, fileMag, filePh)
% Pattern examples at two NMSE levels: two panels, each with the pattern on the left y axis and
% its error on the right y axis. Axis limits are chosen so the pattern sits at the top and the
% error at the bottom, close but not overlapping. Truth solid; the two examples dashed and dotted.
outDir = cpr_paths().results;
E = {};
for d = dir(fullfile(outDir, [study '_w*.mat'])).'
    S = load(fullfile(d.folder, d.name));
    for e = S.res(:).', if ~isempty(e.Xh), E{end+1} = e; end, end %#ok<AGROW>
end
nm = cellfun(@(e) e.raw, E);
[~, i1] = min(abs(nm + 20)); [~, i2] = min(abs(nm + 34)); E = E([i1 i2]);   % the -20 and -34 dB examples
c = E{1}.cfg; [X, ~] = pattern_crlh(c.nw, c.nd, c.f0idx); k = ceil(c.nd/2);
A = X(:, :, k); [~, ip] = max(abs(A(:))); [~, jp] = ind2sub(size(A), ip);
th = (0:c.nw-1)*360/c.nw; a = A(:, jp); pk = max(abs(a));
blue = F.col.IS; red = [0.8 0 0]; ls = {'--', ':'}; mn = char(8722);
leg = {'Truth', sprintf('%s%.0f dB', mn, -E{1}.raw), sprintf('%s%.0f dB', mn, -E{2}.raw), ...
       sprintf('error, %s%.0f dB', mn, -E{1}.raw), sprintf('error, %s%.0f dB', mn, -E{2}.raw)};
B = cellfun(@(e) double(e.Xh(:, jp, k)), E, 'UniformOutput', false);
% ---- magnitude
db = @(v) 20*log10(max(abs(v), 1e-6)/pk);
L = [db(a), db(B{1}), db(B{2})];
Er = [db(abs(B{1}) - abs(a)), db(abs(B{2}) - abs(a))];          % magnitude error level re the peak
Er = max(Er, -80);
twoaxis(F, th, L, Er, blue, red, ls, leg, 'magnitude (dB)', 'magnitude error (dB)', fileMag);
% ---- phase (truth unwrapped along the cut; recovered = truth + phase error)
pT = unwrap(angle(a))*180/pi;
pe = [angle(B{1} .* conj(a)), angle(B{2} .* conj(a))]*180/pi;
L = [pT, pT + pe(:, 1), pT + pe(:, 2)];
twoaxis(F, th, L, pe, blue, red, ls, leg, 'phase (deg)', 'phase error (deg)', filePh);
end

function twoaxis(F, th, L, Er, blue, red, ls, leg, ylL, ylR, file)
f = F.new(); ax = axes(f); hold(ax, 'on'); h = gobjects(1, 5);
yyaxis(ax, 'left');
h(1) = plot(ax, th, L(:, 1), '-', 'Color', 'k', 'LineWidth', F.lw);
h(2) = plot(ax, th, L(:, 2), ls{1}, 'Color', blue, 'LineWidth', F.lw);
h(3) = plot(ax, th, L(:, 3), ls{2}, 'Color', blue, 'LineWidth', F.lw);
lo = min(L(:)); hi = max(L(:)); H = (hi - lo)/0.42;               % pattern in 55-97% of the height
ylim(ax, [hi + 0.03*H - H, hi + 0.03*H]); ylabel(ax, ylL); ax.YAxis(1).Color = blue;
yyaxis(ax, 'right');
h(4) = plot(ax, th, Er(:, 1), ls{1}, 'Color', red, 'LineWidth', F.lw);
h(5) = plot(ax, th, Er(:, 2), ls{2}, 'Color', red, 'LineWidth', F.lw);
lo = min(Er(:)); hi = max(Er(:)); H = (hi - lo)/0.27;             % error in 3-30% of the height
ylim(ax, [lo - 0.03*H, lo - 0.03*H + H]); ylabel(ax, ylR); ax.YAxis(2).Color = red;
xlim(ax, [0 360]); xlabel(ax, '\theta (deg)');
F.style(ax);
% legend in the gap between the two bands (explicit position; 'best' is unreliable here)
lg = legend(ax, h, leg, 'NumColumns', 3, 'FontSize', 7); drawnow;
lp = lg.Position; axp = ax.Position;
lg.Position = [axp(1) + (axp(3) - lp(3))/2, axp(2) + 0.425*axp(4) - lp(4)/2, lp(3), lp(4)];
F.save(f, file);
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
