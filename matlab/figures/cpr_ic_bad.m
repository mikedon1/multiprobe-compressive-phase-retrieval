function r = cpr_ic_bad(seeds, blocks)
%CPR_IC_BAD  Phase retrieval with the unmodified IC apparatus (no leg probes), paper Fig. IC_bad.
%   cpr_ic_bad                tries seeds 1-3 and plots the trial whose magnitude is recovered well and
%                             whose phase differs most between the plotted cuts
%   cpr_ic_bad(seeds, blocks) tries the given seeds and plots the given phi positions (cuts)
%The IC base configuration of the current preset (pattern size, wavelengths, K = 4 tones, sampling
%ratio, restarts) at 100 dB SNR, but with nlegs = 0. Only the first wavelength is plotted. The magnitude
%is recovered but each cut in phi keeps its own phase factor, so the phase is wrong. Four consecutive
%cuts (blocks = 1:4) are plotted, vectorized, with vertical lines between them and the phase
%unwrapped. One global phase is removed from the recovery.
%Writes <results>/figs/paper/IC_bad.pdf/.png when cpr_settings().figures.print is true.
if nargin < 1, seeds = 1:3; end
if nargin < 2, blocks = 1:4; end
c = cpr_cfg('arch', 'IC', 'nlegs', 0, 'snrdB', 100);
[X, ~] = pattern_crlh(c.nw, c.nd, c.f0idx);
x = reshape(X(:, :, 1), [], 1); nb = c.nw;                 % first wavelength
idx = reshape((blocks(:)' - 1)*nb + (1:nb)', [], 1);       % entries of the plotted cuts
best = -inf; sel = 0;
for s = seeds
    rs = cpr_run(c, s); xh = double(reshape(rs.Xh(:, :, 1), [], 1));
    xa = xh * exp(1i*angle(xh' * x));                      % one global phase only
    magErr = 10*log10(sum((abs(xh) - abs(x)).^2)/sum(abs(x).^2));
    off = zeros(1, numel(blocks));                         % phase factor of each plotted cut (deg)
    for k = 1:numel(blocks), i = (blocks(k)-1)*nb + (1:nb); off(k) = angle(x(i)' * xa(i)) * 180/pi; end
    spread = max(abs(off - mean(off)));
    fprintf('seed %d: magnitude error %.1f dB, phase offsets of the plotted cuts [%s] deg\n', s, magErr, num2str(off, '%.0f '));
    if magErr < -30 && spread > best, best = spread; r = rs; xbest = xa; sel = s; end
end
if sel == 0, error('cpr_ic_bad: no seed gave a good magnitude'); end
fprintf('plotting seed %d\n', sel);
xt = x(idx); xr = xbest(idx); n = numel(idx); k = 1:n;
F = cpr_fig(); f = F.new(12, 6);
ax1 = subplot(2, 1, 1); hold(ax1, 'on');
h1 = plot(ax1, k, unwrap(angle(xr)), '-', 'Color', F.col.IS, 'LineWidth', F.lw);
h2 = plot(ax1, k, unwrap(angle(xt)), '--', 'Color', F.col.RS, 'LineWidth', F.lw);
ylabel(ax1, 'Phase'); F.style(ax1);
yl = ylim(ax1); t = ceil(2*yl(1)/pi):floor(2*yl(2)/pi); t = t/2;     % ticks at multiples of pi/2
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
% fixed layout: same plot width in both panels, room for the x label
ax1.Position = [0.11 0.60 0.60 0.36]; ax2.Position = [0.11 0.19 0.60 0.36];
pDir = fullfile(cpr_paths().results, 'figs', 'paper'); if ~exist(pDir, 'dir'), mkdir(pDir); end
F.save(f, fullfile(pDir, 'IC_bad'));
end

function s = piLabel(m)
if m == 0, s = '0'; return, end
q = round(2*m);                                   % m in halves of pi
if mod(q, 2) == 0, k = q/2; if abs(k) == 1, s = [repmat('-', 1, k < 0) '\pi']; else, s = sprintf('%d\\pi', k); end
else, if abs(q) == 1, s = [repmat('-', 1, q < 0) '\pi/2']; else, s = sprintf('%d\\pi/2', q); end
end
end
