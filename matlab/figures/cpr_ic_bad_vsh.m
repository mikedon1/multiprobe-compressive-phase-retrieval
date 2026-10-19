function cpr_ic_bad_vsh(blocks)
%CPR_IC_BAD_VSH  Paper Fig. 3(c) from the saved pr_icbad run (2026-10-19): IC with the 2D Fourier basis and no
%leg probes at the paper's settings (32 x 16, both components, m/n 0.3). Same figure as cpr_ic_bad.m (which runs
%its own recovery with the current defaults, now the spherical-wave basis), but read from the study results
%(<results>/pr_icbad_w*.mat, Xh saved). First wavelength, component along the scan circle, the cuts in phi
%given by blocks (default 1:4), vectorized with vertical lines between them, phase unwrapped, one global phase
%removed. Writes <results>/figs/paper/IC_bad.pdf.
if nargin < 1, blocks = 1:4; end
E = [];
for d = dir(fullfile(cpr_paths().results, 'pr_icbad_w*.mat')).'
    S0 = load(fullfile(d.folder, d.name));
    for e = S0.res(:).', if isfield(e, 'Xh') && ~isempty(e.Xh), E = e; end, end
end
if isempty(E), error('cpr_ic_bad_vsh: no pr_icbad result with Xh'); end
c = E.cfg; [X, ~] = pattern_crlh_vsh(c.nw, c.nd, c.f0idx, 16, c.npol);
x = reshape(X(:, :, 1), [], 1); xh = double(reshape(E.Xh(:, :, 1), [], 1)); nb = c.nw;
xa = xh * exp(1i*angle(xh' * x));
idx = reshape((blocks(:)' - 1)*nb + (1:nb)', [], 1);
off = arrayfun(@(b) angle(x((b-1)*nb + (1:nb))' * xa((b-1)*nb + (1:nb))) * 180/pi, blocks);
fprintf('pr_icbad seed %d: NMSE %.1f dB (global), magnitude error %.1f dB, phase offsets of the cuts [%s] deg\n', ...
    E.seed, E.global, 10*log10(sum((abs(xh) - abs(x)).^2)/sum(abs(x).^2)), num2str(off, '%.0f '));
xt = x(idx); xr = xa(idx); n = numel(idx); k = 1:n;
F = cpr_fig(); f = F.new(12, 6);
ax1 = subplot(2, 1, 1); hold(ax1, 'on');
h1 = plot(ax1, k, unwrap(angle(xr)), '-', 'Color', F.col.IS, 'LineWidth', F.lw);
h2 = plot(ax1, k, unwrap(angle(xt)), '--', 'Color', F.col.RS, 'LineWidth', F.lw);
ylabel(ax1, 'Phase'); F.style(ax1);
yl = ylim(ax1); t = ceil(2*yl(1)/pi):floor(2*yl(2)/pi); t = t/2;
if numel(t) > 6, t = t(mod(t, 1) == 0); end
set(ax1, 'YTick', t*pi, 'YTickLabel', arrayfun(@piLabel, t, 'UniformOutput', false), 'XTickLabel', []);
legend(ax1, [h1 h2], {'$\angle\hat{\mathbf{x}}$', '$\angle\mathbf{x}$'}, 'Interpreter', 'latex', ...
    'Location', 'northeastoutside', 'AutoUpdate', 'off');
ax2 = subplot(2, 1, 2); hold(ax2, 'on');
h3 = plot(ax2, k, abs(xr), '-', 'Color', F.col.IS, 'LineWidth', F.lw);
h4 = plot(ax2, k, abs(xt), '--', 'Color', F.col.RS, 'LineWidth', F.lw);
ylabel(ax2, 'Magnitude'); xlabel(ax2, 'Index'); F.style(ax2);
legend(ax2, [h3 h4], {'$|\hat{\mathbf{x}}|$', '$|\mathbf{x}|$'}, 'Interpreter', 'latex', ...
    'Location', 'northeastoutside', 'AutoUpdate', 'off');
for ax = [ax1 ax2]
    xlim(ax, [1 n]); set(ax, 'TickLabelInterpreter', 'tex');
    for b = 1:numel(blocks)-1, xline(ax, b*nb, '-', 'Color', 'k', 'LineWidth', 0.75); end
end
ax1.Position = [0.11 0.60 0.60 0.36]; ax2.Position = [0.11 0.19 0.60 0.36];
pDir = fullfile(cpr_paths().results, 'figs', 'paper'); if ~exist(pDir, 'dir'), mkdir(pDir); end
F.save(f, fullfile(pDir, 'IC_bad'));
end

function s = piLabel(m)
if m == 0, s = '0'; return, end
q = round(2*m);
if mod(q, 2) == 0, k = q/2; if abs(k) == 1, s = [repmat('-', 1, k < 0) '\pi']; else, s = sprintf('%d\\pi', k); end
else, if abs(q) == 1, s = [repmat('-', 1, q < 0) '\pi/2']; else, s = sprintf('%d\\pi/2', q); end
end
end
