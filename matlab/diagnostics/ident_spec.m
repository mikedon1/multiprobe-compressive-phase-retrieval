function ident_spec()
% Full spectrum at L = 12 of the normalized out-of-span fraction for (a) per-rotation phases, (b) IC per-port
% cable lengths; plus the explicit smooth candidate dL ~ cos(theta_physical), same for both pols.
cd(fileparts(fileparts(fileparts(mfilename('fullpath'))))); setup; cpr_settings('paper');
maxNumCompThreads(1); warning('off', 'all');
c0 = cpr_cfg('arch', 'IS');
[X, f] = pattern_crlh_vsh(c0.nw, c0.nd, c0.f0idx, 16, c0.npol);
[nw, nh, ns] = size(X); nd = ns / c0.npol; npol = c0.npol; n2 = nw*nh;
stack = @(v) reshape(permute(reshape(v, n2, nd, npol), [1 3 2]), [], nd);
Xs = stack(X(:));
r = (1:npol*n2).'; rr = mod(r-1, n2) + 1; th = mod(rr-1, nw) + 1; col = floor((rr-1)/nw) + 1; pol = floor((r-1)/n2) + 1;
lam = 299792458 ./ f(:).'; lamc = 299792458 / mean(f);
Pb = vsh_basis(nw, nh, 12, npol); out = @(V) V - Pb*(Pb'*V);
Mr = zeros(numel(Xs), nh); Mo = Mr;
for l = 1:nh, V = 1i*Xs .* (col == l); Mr(:, l) = V(:); W = out(V); Mo(:, l) = W(:); end
s = spec(Mo, Mr, nh); fprintf('rotation phases, L=12, all %d values (dB): %s\n', numel(s), mat2str(round(20*log10(s)', 1)));
np = nw*npol; Cr = zeros(numel(Xs), np); Co = Cr;
for p = 1:np
    t = mod(p-1, nw) + 1; u = floor((p-1)/nw) + 1;
    V = 1i*(Xs .* ((th == t) & (pol == u))) .* (lamc ./ lam); Cr(:, p) = V(:); W = out(V); Co(:, p) = W(:);
end
s = spec(Co, Cr, np); sd = 20*log10(s');
fprintf('IC cable lengths, L=12, all %d values (dB): %s\n', numel(s), mat2str(round(sd, 0)));
fprintf('  count below -60 dB: %d, below -30 dB: %d, below -10 dB: %d\n', sum(sd < -60), sum(sd < -30), sum(sd < -10));
% explicit smooth candidates: dL_p = cos(k * theta_phys(p)), same for both pols, and theta-hat only
thS = (0:nw-1)' * 360/nw; thP = min(thS, 360 - thS);
for kk = 1:3
    g = repmat(cosd(kk*thP), npol, 1); g = g - mean(g);
    fprintf('  dL ~ cos(%d theta), both pols: %.1f dB\n', kk, 20*log10(norm(Co*g)/norm(Cr*g)));
end
g = [cosd(thP); zeros(nw,1)]; g = g - mean(g); fprintf('  dL ~ cos(theta), theta-hat ports only: %.1f dB\n', 20*log10(norm(Co*g)/norm(Cr*g)));
g = [ones(nw,1); -ones(nw,1)]; fprintf('  +d on along ports, -d on across ports: %.1f dB\n', 20*log10(norm(Co*g)/norm(Cr*g)));
disp('SPEC DONE')
end
function s = spec(Mout, Mraw, k)
Q = null(ones(1, k)); Wo = [real(Mout); imag(Mout)] * Q; Wr = [real(Mraw); imag(Mraw)] * Q;
[~, Rw] = qr(Wr, 0); s = sort(svd(Wo / Rw));
end
