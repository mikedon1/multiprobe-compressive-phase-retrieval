function ident_margins()
% Identifiability margins of the VSH span (no new simulation; linear algebra on the example truth).
% For L = 8..20 on the paper grid (32 x 16, dual pol, 16 wavelengths):
%  rnd  : energy outside the span of x with a random phase per rotation (phi column), 3 draws (as t7)
%  rotW : worst case for SMALL per-rotation phases g (g perp. constant): min ||(I-P) i x.*g|| / ||x|| per
%         1 rad rms of g, i.e. smallest nontrivial singular value of the linearized map, in dB
%  rot2 : second smallest (to see whether one direction or many are weak)
%  cabW : same for IC per-port cable LENGTH errors (64 ports = 32 ring positions x 2 pols), phase
%         2*pi*dL/lambda_k, dL perp. constant, per 1 deg rms phase at band centre
cd(fileparts(fileparts(fileparts(mfilename('fullpath'))))); setup; cpr_settings('paper');
maxNumCompThreads(1); warning('off', 'all');
c0 = cpr_cfg('arch', 'IS');
[X, f] = pattern_crlh_vsh(c0.nw, c0.nd, c0.f0idx, 16, c0.npol);
[nw, nh, ns] = size(X); nd = ns / c0.npol; npol = c0.npol; n2 = nw*nh;
stack = @(v) reshape(permute(reshape(v, n2, nd, npol), [1 3 2]), [], nd);   % (npol*n2) x nd
Xs = stack(X(:)); nx = norm(Xs, 'fro');
r = (1:npol*n2).'; rr = mod(r-1, n2) + 1; th = mod(rr-1, nw) + 1; col = floor((rr-1)/nw) + 1; pol = floor((r-1)/n2) + 1;
lam = 299792458 ./ f(:).'; lamc = 299792458 / mean(f);
fprintf('nw %d nh %d npol %d nd %d\n', nw, nh, npol, nd);
fprintf('%4s %5s %18s %9s %9s %9s %9s\n', 'L', 'cols', 'rnd outside (dB)', 'rotW', 'rot2', 'cabW', 'cab2');
for L = [8 10 12 14 16 18 20]
    Pb = vsh_basis(nw, nh, L, npol);
    out = @(V) V - Pb*(Pb'*V);
    rng(1); R = zeros(1,3);
    for s = 1:3
        ph = exp(1i*2*pi*rand(1, nh)); V = Xs .* ph(col).'; R(s) = 20*log10(norm(out(V), 'fro')/norm(V, 'fro'));
    end
    % rotation phases: column l of M = (I-P)(i x masked to column l), all wavelengths
    M = zeros(numel(Xs), nh);
    for l = 1:nh, V = 1i*Xs .* (col == l); M(:, l) = reshape(out(V), [], 1); end
    Mraw = zeros(numel(Xs), nh); for l = 1:nh, Mraw(:, l) = reshape(1i*Xs .* (col == l), [], 1); end; rotSV = gmargin(M, Mraw, nh);
    % IC cable lengths: port = (theta row, pol); phase 2*pi*dL/lambda_k; dL in units giving 1 deg rms at centre
    np = nw*npol; C = zeros(numel(Xs), np);
    for p = 1:np
        t = mod(p-1, nw) + 1; u = floor((p-1)/nw) + 1; mask = (th == t) & (pol == u);
        V = 1i*(Xs .* mask) .* (lamc ./ lam) * (pi/180);       % d(phase) = (lamc/lam_k) * 1 deg per unit dL
        C(:, p) = reshape(out(V), [], 1);
    end
    Craw = zeros(numel(Xs), np); for p = 1:np, t = mod(p-1, nw) + 1; u = floor((p-1)/nw) + 1; Craw(:, p) = reshape(1i*(Xs .* ((th == t) & (pol == u))) .* (lamc ./ lam), [], 1); end; e = sum(abs(Craw).^2, 1); [~, iw] = sort(e); if L == 12, fprintf("weakest ports (row,pol): %s  energy share %s\n", mat2str([mod(iw(1:4)-1,nw)+1; floor((iw(1:4)-1)/nw)+1]), mat2str(e(iw(1:4))/sum(e), 2)); end; cabSV = gmargin(C, Craw, np);
    fprintf('%4d %5d %18s %9.1f %9.1f %9.1f %9.1f\n', L, size(Pb,2), mat2str(round(R,1)), ...
        20*log10(rotSV(1)), 20*log10(rotSV(2)), 20*log10(cabSV(1)), 20*log10(cabSV(2)));
end
disp('IDENT DONE')
end

function sv = margin(M, k)
% smallest two singular values of the real-linear map g -> M g restricted to g perp. ones, per unit RMS of g
Q = null(ones(1, k));                         % k x (k-1), orthonormal, perp. to constant
Mr = [real(M); imag(M)] * Q;
s = svd(Mr, 'econ'); s = sort(s);
sv = s(1:2) * sqrt(k);                         % unit-norm g has rms 1/sqrt(k): rescale to unit rms
end

function sv = gmargin(Mout, Mraw, k)
% smallest two values of ||Mout g|| / ||Mraw g|| over real g perp. ones: the fraction of a perturbation's
% own amplitude that lies outside the span (0 dB = all outside, very negative = invisible to the basis)
Q = null(ones(1, k));
Wo = [real(Mout); imag(Mout)] * Q; Wr = [real(Mraw); imag(Mraw)] * Q;
[~, Rw] = qr(Wr, 0);
s = sort(svd(Wo / Rw)); sv = s(1:2);
end
