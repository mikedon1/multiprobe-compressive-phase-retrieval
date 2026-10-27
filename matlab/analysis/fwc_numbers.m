function fwc_numbers()
%FWC_NUMBERS  Numbers for the Section VII text and Table S6 from the finished HFSS ring-calibration trials.
cpr_setup(); maxNumCompThreads(1); warning('off', 'all');
V = load(fullfile(cpr_paths().data, 'fullwave', 'fwc_validation.mat')); V = V.V; a = V.A.a; f = V.freq(:).'; nf = numel(f);
Tr = fwc_ring_data('hfss');
[T, C] = fwc_load_trials('^fwc_trials_true_hfss_ring\.mat$');

ref = arrayfun(@(t) t.nmse.perfreqScale, T); raw = arrayfun(@(t) t.nmse.rawPerfreqScale, T);
fprintf('TABLE S6 ring columns (seed, PR-GAMP, refined):\n'); for k = 1:numel(T), fprintf('  %d & %.2f & %.2f\n', T(k).seed, raw(k), ref(k)); end
lin = @(x) 10*log10(mean(10.^(x/10)));
fprintf('  mean of dB: %.2f %.2f | linear: %.2f %.2f | refined best %.2f worst %.2f\n', mean(raw), mean(ref), lin(raw), lin(ref), min(ref), max(ref));
% factor deviations from a constant and one delay, benchmark under the same alignment
ampdev = zeros(1, numel(T)); phres = ampdev; bench = ampdev; gd = ampdev; tres = ampdev;
for k = 1:numel(T)
    kf = T(k).perfreqFactor(:).'; ph = unwrap(angle(kf)); p = polyfit(f - f(1), ph, 1); res = ph - polyval(p, f - f(1));
    phres(k) = rms(res) * 180/pi; ampdev(k) = 20*log10(rms(abs(kf)/mean(abs(kf)) - 1) + 1) ; ampdev(k) = rms(20*log10(abs(kf)/mean(abs(kf))));
    ahRef = a ./ kf;                                         % reference coefficients before the fitted factors
    ahD = ahRef .* exp(1i * polyval(p, f - f(1)));
    bench(k) = scoreScale(a(:), ahD(:));
    gd(k) = T(k).nmse.globalDelay;
end
fprintf('FACTOR: amplitude rms %.2f dB, phase residual after constant+delay %.2f deg (means over trials)\n', mean(ampdev), mean(phres));
fprintf('BENCHMARK constant+delay: recovered mean %.2f (best %.2f worst %.2f) | reference benchmark mean %.2f | difference %.2f dB\n', mean(gd), min(gd), max(gd), mean(bench), mean(gd) - mean(bench));
% residual of the true coefficients against the noisy amplitudes, same readings as the trials
D = struct('Xring', Tr.Xring, 'R', Tr.Rring, 'modes', V.A.modes, 'psi', Tr.psi, 'f', f, 'aTrue', a);
[np, nrot, ~] = size(D.Xring); nm = size(D.modes, 1); G = cell(1, nf);
for q = 1:nf
    E = exp(-1i * D.modes(:, 2).' .* D.psi(:)); Gf = zeros(np*nrot, nm);
    for l = 1:nrot, Gf((l-1)*np + (1:np), :) = D.R(:, :, q) .* E(l, :); end
    G{q} = Gf;
end
for k = 1:numel(T)
    c = C; c.nlegs = 0; c = copyOver(c, T(k));
    Dm = struct('Xring', D.Xring, 'f', f, 'nlegsAvail', 0, 'Xleg', [], 'anchor', 1, 'legShift', 1);
    Sd = make_readings_fw(c, Dm, T(k).seed, 'ideal'); rng(T(k).seed + 0);                       % same readings as fwc_recover
    z = Sd.Phi * D.Xring(:); wvar = cpr_noisevar(c, z, D.Xring);
    Xm = zeros(np*nrot, nf); for q = 1:nf, Xm(:, q) = G{q} * a(:, q); end
    Xd = reshape(D.Xring, np*nrot, nf); for q = 1:nf, Xm(:, q) = Xm(:, q) * ((Xm(:, q)' * Xd(:, q)) / (Xm(:, q)' * Xm(:, q))); end
    zm = Sd.Phi * Xm(:);
    tres(k) = 20*log10(norm(abs(zm) - abs(z)) / norm(abs(z)));
    nz = 10*log10(wvar / mean(abs(z).^2));
    if k == 1, fprintf('noise level relative to the signal: %.1f dB\n', nz); end
end
fprintf('TRUTH amplitude residual vs the noise-free HFSS readings (model error only): mean %.2f dB (%s); accepted PR-GAMP residual mean %.2f dB\n', mean(tres), mat2str(round(tres,1)), mean(arrayfun(@(t) t.nres_dB, T)));
% other sets
names = {'fwc_trials_true_hfss_ring_r60','fwc_trials_true_az0.5_rel-Inf_hfss_ring_pa05','fwc_trials_true_az1_rel-Inf_hfss_ring_pa1','fwc_trials_true_az0_rel-40_hfss_ring_pr40','fwc_trials_true_az0_rel-30_hfss_ring_pr30'};
lab = {'m/n 0.6','az 0.5 deg','az 1 deg','entries -40 dB','entries -30 dB'};
for i = 1:numel(names)
    S = load(fullfile(cpr_paths().results, [names{i} '.mat'])); v = arrayfun(@(t) t.nmse.perfreqScale, S.T); r = arrayfun(@(t) t.nmse.rawPerfreqScale, S.T);
    fprintf('SET %-15s seeds %s: refined %s mean %.2f | PR-GAMP alone mean %.2f | constant+delay mean %.2f\n', lab{i}, mat2str([S.T.seed]), mat2str(round(v,2)), mean(v), mean(r), mean(arrayfun(@(t) t.nmse.globalDelay, S.T)));
end
fprintf('unperturbed seeds 1-3: refined %.2f, PR-GAMP %.2f\n', mean(ref(1:3)), mean(raw(1:3)));
end

function c = copyOver(c, t) %#ok<INUSD>
end

function e = scoreScale(t, y)
k = (y' * t) / max(y' * y, eps); e = 10*log10(norm(k*y - t)^2 / norm(t)^2);
end