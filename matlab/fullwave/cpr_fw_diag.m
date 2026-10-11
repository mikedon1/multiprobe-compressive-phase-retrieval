function cpr_fw_diag()
%CPR_FW_DIAG  Section VI data checks on the 32-rotation multi-frequency HFSS sweep (data/hfss):
%per-leg mismatch against the anchor mapping, the 99%-energy coefficient count of the ring pattern,
%and the frequency with the most detailed pattern. Saves <results>/fw_diag.mat and
%fw_data_hfss.mat (the data structure D used by the trial runs).
outDir = cpr_paths().results;
D = fw_data('hfss', 'Glob', fullfile(cpr_paths().data, 'hfss', 'sweep_redesign_mf_rot*.s20p'));
save(fullfile(outDir, 'fw_data_hfss.mat'), 'D');
[nTh, nPh, nf] = size(D.Xring);
fprintf('%s: %d x %d x %d, f = %.3f..%.3f GHz, legShift %+d (residual %.3f, other direction %.3f)\n', ...
    D.src, nTh, nPh, nf, D.f(1)/1e9, D.f(end)/1e9, D.legShift, D.legResidual, D.legResidualOther);
% per-leg mismatch: leg k at rotation L vs anchor at L + legShift*k, all frequencies
a = squeeze(D.Xring(D.anchor, :, :));                       % nPh x nf
mm = zeros(1, D.nlegsAvail);
for k = 1:D.nlegsAvail
    l = squeeze(D.Xleg(k, :, :));
    ref = circshift(a, -D.legShift*k, 1);                   % ref(L) = a(L + legShift*k)
    mm(k) = norm(l - ref, 'fro') / norm(l, 'fro');
    fprintf('leg %d: relative difference from the anchor mapping %.1f%% (%.1f dB)\n', k, 100*mm(k), 20*log10(mm(k)));
end
% 99%-energy coefficient count, joint 2D DFT x delay-compensated DCT (Eq. (9) delay)
X2 = reshape(fft2(D.Xring) / sqrt(nTh*nPh), nTh*nPh, nf);
step = angle(sum(sum(conj(X2(:, 1:end-1)) .* X2(:, 2:end))));
tau = -step / (2*pi*(D.f(2) - D.f(1)));
S2 = X2 * conj(cpr_basis('dct', D.f, tau));
e = sort(abs(S2(:)).^2, 'descend'); K99 = find(cumsum(e) >= 0.99*sum(e), 1);
S0 = X2 * conj(cpr_basis('dct', D.f, 0));
e0 = sort(abs(S0(:)).^2, 'descend'); K99nodelay = find(cumsum(e0) >= 0.99*sum(e0), 1);
fprintf('K99 = %d of %d (delay %.2f ns); without delay compensation %d\n', K99, numel(D.Xring), 1e9*tau, K99nodelay);
% "most detailed" frequency: most 2D DFT coefficients for 99% at that frequency alone
kf = zeros(1, nf);
for t = 1:nf
    c = sort(abs(X2(:, t)).^2, 'descend'); kf(t) = find(cumsum(c) >= 0.99*sum(c), 1);
end
[~, tmax] = max(kf);
fprintf('per-frequency 99%% counts: %s -> most detailed %.3f GHz\n', mat2str(kf), D.f(tmax)/1e9);
save(fullfile(outDir, 'fw_diag.mat'), 'mm', 'K99', 'K99nodelay', 'tau', 'kf', 'tmax');
end
