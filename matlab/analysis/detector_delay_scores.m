function detector_delay_scores(part)
%DETECTOR_DELAY_SCORES  Re-runs the low-pass points of pr_detbw and scores them with one global phase, one phase per wavelength AND one
%constant plus one delay for the whole band (a pure group delay of the detector is indistinguishable from the AUT's own delay).
%   part 1 or 2 splits the 20 runs between two workers. Output: results_extra/detdelay_p<part>.mat and printed lines.
cpr_setup(); cpr_settings('paper'); maxNumCompThreads(2);

cfgs = {};
for a = {'IS', 'IC'}, for sp = [16 8], for bw = [4 2 1 0.5 0.25]
    cfgs(end+1, :) = {sprintf('%s_span%d_bw%g', a{1}, sp, bw), {'arch', a{1}, 'toneSpan', sp, 'detBwGHz', bw, 'maxTry', 3}, bw}; %#ok<AGROW>
end, end, end
% ideal references for the same settings are in pr_detbw (bwInf); add them here too with the same scores
for a = {'IS', 'IC'}, for sp = [16 8]
    cfgs(end+1, :) = {sprintf('%s_span%d_bwInf', a{1}, sp), {'arch', a{1}, 'toneSpan', sp, 'maxTry', 3}, Inf}; %#ok<AGROW>
end, end
idx = part:2:size(cfgs, 1); out = struct('label', {}, 'global', {}, 'perlambda', {}, 'delay', {}, 'tauPs', {}, 'tauTheoryPs', {});
for i = idx
    c = cpr_cfg(cfgs{i, 2}{:}); r = cpr_run(c, 1);
    [X, f] = pattern_crlh_vsh(c.nw, c.nd, c.f0idx, 16, c.npol); nd = c.nd; Xh = double(reshape(r.Xh, size(X)));
    kf = zeros(1, nd);
    for k = 1:nd, cols = [k, nd + k]; y = reshape(Xh(:, :, cols), [], 1); t = reshape(X(:, :, cols), [], 1); kf(k) = (y' * t) / max(y' * y, eps); end
    ph = unwrap(angle(kf)); p = polyfit(f(:).' - f(1), ph, 1);
    Xa = Xh; for k = 1:nd, Xa(:, :, [k nd+k]) = Xh(:, :, [k nd+k]) * exp(1i * polyval(p, f(k) - f(1))); end
    y = Xa(:); t = X(:); kc = (y' * t) / (y' * y); de = 20*log10(norm(kc*y - t) / norm(t));
    e = struct('label', cfgs{i, 1}, 'global', r.score.global, 'perlambda', r.score.perlambda, 'delay', de, 'tauPs', p(1)/(2*pi)*1e12, ...
        'tauTheoryPs', 1e12/(2*pi*cfgs{i, 3}*1e9));
    out(end+1) = e; %#ok<AGROW>
    fprintf('%-20s global %7.2f | per-wavelength %7.2f | constant+delay %7.2f dB | fitted delay %7.1f ps (low-pass group delay %5.1f ps)\n', ...
        e.label, e.global, e.perlambda, e.delay, e.tauPs, e.tauTheoryPs);
    save(fullfile(cpr_paths().results, sprintf('detdelay_p%d.mat', part)), 'out');
end
end