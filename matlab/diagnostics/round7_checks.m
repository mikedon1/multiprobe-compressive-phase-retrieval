function round7_checks(part)
%ROUND7_CHECKS  Short checks for the round-7 critique (2026-10-20). part = 1..5 (default all).
% 1  T1 counterexample E+- = -sin^2(theta) e^{+-i phi} (0.4 theta_hat + phi_hat) on the paper scan grid: in the VSH
%    span at L = 2 and 12? each cut related by e^{-2i phi}?
% 2  structured rotation phases g_l = c*phi_l (c = +-1..+-6) applied to the example truth (CRLH, and the 5x5 array on
%    64 x 32 with L = 20): fraction outside the span (finite phases, not small)
% 3  C4: Gram-Schmidt R between modes of different m on the scan grid (does e^{-i m psi} hold in Q coordinates?)
% 4  T2: compressed-regime local tests at the base readings (IC, m/n 0.3, K 4, 16 freqs), coefficient support =
%    largest per-frequency VSH coefficients holding 99.9% / 99.99% of the energy:
%    (a) amplitude Jacobian on the support: rank deficiency beyond the global phase
%    (b) rotation phases: fraction of i g.x outside span(support columns), worst g (g perp 1)
%    (c) cable lengths: (I - P_Js) J_L, smallest normalized singular values
% 5  actual m for IS/RS/IC/RC at the base values
if nargin < 1, part = 1:5; end
cd(fileparts(fileparts(fileparts(mfilename('fullpath'))))); setup; cpr_settings('paper');
maxNumCompThreads(1); warning('off', 'all');
c0 = cpr_cfg(); nw = c0.nw; nh = c0.nh; n2 = nw*nh;
[th, phs, sg] = vsh_scan(nw, nh);              % true directions, source-frame azimuth, far-half sign
ts = mod((0:n2-1)', nw) * 2*pi/nw; far = ts > pi + 1e-9;
phTrue = mod(phs + pi/2, 2*pi);                % pattern_crlh rotates the source azimuth by 90 deg
col = floor((0:n2-1)'/nw) + 1;

if any(part == 1)
    for s = [1 -1]
        Et = -sin(th).^2 .* exp(1i*s*phTrue) * 0.4; Ep = -sin(th).^2 .* exp(1i*s*phTrue);
        % scan components: along = sg*theta, across = sg*phi (vsh_basis npol = 2 order [along; across])
        X{(3-s)/2} = [sg .* Et; sg .* Ep]; %#ok<AGROW>
    end
    for L = [2 12]
        P = vsh_basis(nw, nh, L, 2);
        r = cellfun(@(v) norm(v - P*(P'*v))/norm(v), X);
        fprintf('(1) L = %2d: E+ outside span %.1e, E- outside span %.1e\n', L, r(1), r(2));
    end
    ratio = X{2} ./ X{1}; ok = abs(X{1}) > 1e-6*max(abs(X{1}));
    phcol = repmat(phs(:), 2, 1);
    % per cut (scan column), the ratio must be one constant
    cc = repmat(col, 2, 1); dev = 0;
    for l = 1:nh, rl = ratio(ok & cc == l); dev = max(dev, max(abs(rl - rl(1)))); end
    fprintf('    E-/E+ constant within every cut to %.1e (|ratio| = %.3f)\n', dev, abs(ratio(find(ok, 1))));
end

if any(part == 2)
    for pat = {'crlh', 'array5'}
        if strcmp(pat{1}, 'crlh'), g = nw; L = 12; Lf = 16; else, g = 64; L = 20; Lf = 30; end
        [X, ~] = pattern_crlh_vsh(g, 16, 1, Lf, 2, pat{1});
        [gw, gh, ns] = size(X); nd = ns/2; m2 = gw*gh;
        Xs = reshape(permute(reshape(X, m2, nd, 2), [1 3 2]), 2*m2, nd);
        P = vsh_basis(gw, gh, L, 2); cl = repmat(floor((0:m2-1)'/gw), 2, 1);
        out = zeros(1, 12); cs = [-6:-1 1:6];
        for i = 1:12
            V = Xs .* exp(1i*cs(i)*cl*pi/gh); out(i) = 20*log10(norm(V - P*(P'*V), 'fro')/norm(V, 'fro'));
        end
        fprintf('(2) %s (%dx%d, L %d): ramp c = %s -> outside span %s dB\n', pat{1}, gw, gh, L, mat2str(cs), mat2str(round(out, 1)));
    end
end

if any(part == 3)
    L = 12; [Yt, Yp, modes] = vsh_eval(L, th, phs); Y = [sg .* Yt; sg .* Yp];
    [Q, R] = qr(Y, 0);                          % same nested subspaces as degree-ordered Gram-Schmidt
    mm = modes(:, 2); D = abs(R) .* (mm ~= mm.');
    fprintf('(3) |R| between different m: max %.1e (relative to max|R| %.1e)\n', max(D(:)), max(abs(R(:))));
    Pb = vsh_basis(nw, nh, L, 2); fprintf('    vsh_basis vs QR span difference %.1e\n', norm(Pb - Q*(Q'*Pb), 'fro'));
end

if any(part == 4)
    c = cpr_cfg('arch', 'IC'); [X, f] = pattern_crlh_vsh(c.nw, c.nd, c.f0idx, 16, 2);
    [nw, nh, ns] = size(X); nd = ns/2; n2 = nw*nh; x = X(:);
    S = make_readings(c, f, 1); z = S.Phi * x; u = z ./ abs(z);
    Pb = vsh_basis(nw, nh, c.vshL, 2); p = size(Pb, 2);
    % pattern-vector basis, one VSH block per frequency (column order a + (k-1) n2 + (q-1) nd n2)
    Psi = sparse(n2*nd*2, p*nd);
    for k = 1:nd, for q = 1:2
        rows = (q-1)*nd*n2 + (k-1)*n2 + (1:n2);
        Psi(rows, (k-1)*p + (1:p)) = Pb((q-1)*n2 + (1:n2), :); %#ok<SPRIX>
    end, end
    s = Psi' * x;
    for tol = [0.999 0.9999]
        [~, o] = sort(abs(s).^2, 'descend'); kk = find(cumsum(abs(s(o)).^2) >= tol*norm(s)^2, 1); Ssup = o(1:kk);
        H = S.Phi * Psi(:, Ssup); uH = full(conj(u) .* H); Js = [real(uH), -imag(uH)];
        sv = svd(full(Js)); nul = sum(sv < 1e-9*sv(1));
        % (b) rotation phases against the support span
        B = Psi(:, Ssup); [Qb, ~] = qr(full(B), 0);
        cl = floor(mod((0:numel(x)-1)', n2)/nw) + 1;
        Mr = zeros(numel(x), nh); Mo = Mr;
        for l = 1:nh, v = 1i*x .* (cl == l); Mr(:, l) = v; Mo(:, l) = v - Qb*(Qb'*v); end
        Qg = null(ones(1, nh)); Wr = [real(Mr); imag(Mr)]*Qg; Wo = [real(Mo); imag(Mo)]*Qg;
        [~, Rw] = qr(Wr, 0); sr = sort(svd(Wo / Rw));
        % (c) cable lengths in the data domain, after removing what the support coefficients can explain
        e = S.ent; lamC = 299792458/mean(f)*100;
        dz = 1i*(lamC ./ reshape(S.lam(e.tone), [], 1))*(pi/180) .* e.v .* x(e.j);
        JL = full(sparse(e.i, e.port, real(conj(u(e.i)) .* dz), S.m, e.nPorts));
        [Us, Ss, ~] = svd(Js, 'econ'); ss = diag(Ss); Us = Us(:, ss > ss(1)*1e-10);
        QL = null(ones(1, e.nPorts)); W = JL*QL; Wo2 = W - Us*(Us'*W); [~, RL] = qr(W, 0);
        sl = sort(20*log10(svd(Wo2 / RL)));
        fprintf('(4) IC base, support %.2f%% (%d complex of %d, m = %d): Js null dim %d (1 = global phase only); rotation worst %.1f dB; cable smallest %s dB, count < -30 dB: %d of %d\n', ...
            100*tol, kk, numel(s), S.m, nul, 20*log10(sr(1)), mat2str(round(sl(1:4)', 1)), sum(sl < -30), numel(sl));
    end
end

if any(part == 5)
    for a = {'IS', 'RS', 'IC', 'RC'}
        c = cpr_cfg('arch', a{1}); [~, f] = pattern_crlh_vsh(c.nw, c.nd, c.f0idx, 16, 2);
        S = make_readings(c, f, 1); fprintf('(5) %s: m = %d, readings R = %d, n = %d, m/n = %.4f\n', a{1}, S.m, S.R, S.n, S.m/S.n);
    end
end
disp('ROUND7 CHECKS DONE')
end
