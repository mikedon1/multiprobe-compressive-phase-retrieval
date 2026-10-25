function pole_check()
%POLE_CHECK  (2026-10-19) Is the leg-free result due to the pole samples (seen by every rotation, like the top probe of
%Tena Sanchez 2021) or to the band limit generally? Rotation-phase margin and random-phase out-of-span energy on the
%paper grid (32 x 16, both components, 16 wavelengths, L = 12) with all samples vs with the pole rows (scan theta = 0
%and 180 deg) removed, and with theta within 15 deg of the poles removed (the HFSS ring's coverage).
cd(fileparts(fileparts(fileparts(mfilename('fullpath'))))); setup; cpr_settings('paper');
maxNumCompThreads(1); warning('off', 'all');
c = cpr_cfg(); [X, ~] = pattern_crlh_vsh(c.nw, c.nd, c.f0idx, 16, 2);
[nw, nh, ns] = size(X); nd = ns/2; n2 = nw*nh;
Xs = reshape(permute(reshape(X, n2, nd, 2), [1 3 2]), 2*n2, nd);
Pb = vsh_basis(nw, nh, 12, 2);
r = (1:2*n2).'; rr = mod(r-1, n2) + 1; th = mod(rr-1, nw); col = floor((rr-1)/nw) + 1;   % th = scan index 0..nw-1
thDeg = th*360/nw; polar = min(thDeg, 360 - thDeg);
cases = {'all samples', true(size(r)); 'poles removed', polar > 1e-9; 'within 15 deg of poles removed', polar >= 15 - 1e-9};
for k = 1:3
    keep = cases{k, 2}; B = orth(Pb(keep, :)); Xk = Xs(keep, :); ck = col(keep);
    out = @(V) V - B*(B'*V);
    Mo = zeros(numel(Xk), nh); Mr = Mo;
    for l = 1:nh, V = 1i*Xk .* (ck == l); Mr(:, l) = V(:); W = out(V); Mo(:, l) = W(:); end
    Q = null(ones(1, nh)); Wo = [real(Mo); imag(Mo)]*Q; Wr = [real(Mr); imag(Mr)]*Q;
    [~, Rw] = qr(Wr, 0); s = sort(svd(Wo / Rw));
    rng(1); R = zeros(1, 3);
    for t = 1:3, ph = exp(1i*2*pi*rand(1, nh)); V = Xk .* ph(ck).'; R(t) = 20*log10(norm(out(V), 'fro')/norm(V, 'fro')); end
    fprintf('%-32s rows %4d of %4d, rank %d: small-phase margin %.1f dB (2nd %.1f), random phases outside %s dB\n', ...
        cases{k, 1}, nnz(keep), numel(keep), size(B, 2), 20*log10(s(1)), 20*log10(s(2)), mat2str(round(R, 1)));
end
disp('POLE DONE')
end
