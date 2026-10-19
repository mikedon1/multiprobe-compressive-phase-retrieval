function [Et1, Ep1, info] = twoscan_solve(A1t, A1p, A2t, A2p, S, k, r1, r2, varargin)
%TWOSCAN_SOLVE  Two-scan technique (two spheres), one wavelength (Section VII).
%
%   [Et1, Ep1, info] = twoscan_solve(|Et(r1)|, |Ep(r1)|, |Et(r2)|, |Ep(r2)|, S, k, r1, r2, ...)
%
% Measured magnitudes of both tangential components (dual-polarized probes, every probe read
% individually) on two concentric spheres r1 < r2, on the theta x phi tensor grid of
% S = swe_setup(...). Returns the complex field on sphere r1.
%
% The unknown is the vector of spherical-wave coefficients c (at sphere r1) for degrees
% n <= Nc. Keeping n <= Nc is the spherical counterpart of the aperture constraint (an antenna
% inside a sphere of radius r0 radiates negligible power in degrees n >> k r0), and it is
% REQUIRED for stability (the quasi-Newton variant: with all degrees kept the iteration diverges
% even from the true field).
%
% Method 'joint' (default; error reduction with a joint update): each iteration
%   1. synthesize both spheres from c, keep their phases, impose the measured magnitudes;
%   2. fit both resulting fields and take the weighted least-squares c that best matches both
%      at once (per coefficient: (w1 f1 + w2 conj(g) f2) / (w1 + w2 |g|^2), g = radial ratio).
% Nothing is divided by a small radial factor, so the update is stable and the weighted data
% residual never increases. w_s = 1/||A_s||^2 balances the two spheres.
% Method 'alternate': the literal planar two-plane loop (sphere 1 -> 2 -> back to 1 by 1/g).
% Kept for comparison: unstable on spheres (the inward step amplifies degrees n > k r1).
%
% Options: 'Nc' (default S.N), 'Method' 'joint', 'MaxIter' 1000, 'Tol' 1e-4 and 'Span' 20
% (stop when the data residual improves by less than Tol, relative, over Span iterations),
% 'Restarts' 1 (extra starts use random phases; lowest residual wins), 'Seed' 1, 'InitPhase'
% {phase_t, phase_p} for the first start (tests only). Stopping uses only the measured data.
p = inputParser;
p.addParameter('Nc', S.N); p.addParameter('MaxIter', 1000);
p.addParameter('Method', 'joint');                                % 'joint' | 'alternate' | 'bfgs'
p.addParameter('Tol', 1e-4); p.addParameter('Span', 20); p.addParameter('Restarts', 1);
p.addParameter('Seed', 1); p.addParameter('InitPhase', []);
p.addParameter('InitCoef', []);           % coefficient vector (S indexing) for the first start ('bfgs')
p.addParameter('Polish', false);          % 'bfgs': Levenberg-Marquardt polish of the best start
p.addParameter('StopRes_dB', -inf);       % 'bfgs': stop starting new starts once the residual is below this (dB)
p.addParameter('Optimizer', 'auto');      % 'bfgs': 'fminunc' (Optimization Toolbox) | 'own' (bfgs_min below, no
                                          %   toolbox) | 'auto' (fminunc if the toolbox is available, else own)
p.parse(varargin{:}); o = p.Results;
t0 = tic;
keep = S.n <= o.Nc;
g = swe_radial(S, k, r2) ./ swe_radial(S, k, r1); g(~keep) = 0;          % sphere 1 -> 2
w1 = 1/(norm(A1t(:))^2 + norm(A1p(:))^2); w2 = 1/(norm(A2t(:))^2 + norm(A2p(:))^2);
den = w1 + w2*abs(g).^2; den(~keep) = 1;
gi = zeros(size(g)); gi(keep) = 1 ./ g(keep);
rs = RandStream('mt19937ar', 'Seed', o.Seed);
best = struct('res', inf);
if strcmp(o.Method, 'bfgs')
    [best, nst] = bfgs_solve(A1t, A1p, A2t, A2p, S, g, keep, w1, w2, o, rs);
    [Et1, Ep1] = swe_synth(S, best.c);
    info.iters = best.iters; info.res_dB = 10*log10(best.res); info.hist = best.hist;
    info.restart = best.restart; info.time = toc(t0); info.opts = o; info.c = best.c; info.starts = nst;
    return
end
for t = 1:o.Restarts
    if t == 1 && ~isempty(o.InitPhase), pt = o.InitPhase{1}; pp = o.InitPhase{2};
    elseif t == 1, pt = zeros(size(A1t)); pp = pt;
    else, pt = 2*pi*rand(rs, size(A1t)); pp = 2*pi*rand(rs, size(A1p)); end
    c = swe_fit(S, A1t .* exp(1j*pt), A1p .* exp(1j*pp)); c(~keep) = 0;
    hist = zeros(1, o.MaxIter); it = 0;
    while it < o.MaxIter
        it = it + 1;
        [E1t, E1p] = swe_synth(S, c);
        [E2t, E2p] = swe_synth(S, c .* g);
        r1res = norm(abs([E1t(:); E1p(:)]) - [A1t(:); A1p(:)])^2;
        r2res = norm(abs([E2t(:); E2p(:)]) - [A2t(:); A2p(:)])^2;
        hist(it) = (w1*r1res + w2*r2res) / 2;
        Y1t = A1t .* exp(1j*angle(E1t)); Y1p = A1p .* exp(1j*angle(E1p));
        Y2t = A2t .* exp(1j*angle(E2t)); Y2p = A2p .* exp(1j*angle(E2p));
        switch o.Method
            case 'joint'
                c = (w1*swe_fit(S, Y1t, Y1p) + w2*conj(g).*swe_fit(S, Y2t, Y2p)) ./ den;
            case 'alternate'                 % 1 -> 2 (impose) -> back to 1 by 1/g (impose)
                [F1t, F1p] = swe_synth(S, swe_fit(S, Y2t, Y2p) .* gi);
                c = swe_fit(S, A1t .* exp(1j*angle(F1t)), A1p .* exp(1j*angle(F1p)));
        end
        c(~keep) = 0;
        if it > o.Span && hist(it - o.Span) - hist(it) < o.Tol * hist(it - o.Span), break; end
    end
    if hist(it) < best.res
        best = struct('res', hist(it), 'c', c, 'iters', it, 'hist', hist(1:it), 'restart', t);
    end
end
[Et1, Ep1] = swe_synth(S, best.c);
info.iters = best.iters; info.res_dB = 10*log10(best.res); info.hist = best.hist;
info.restart = best.restart; info.time = toc(t0); info.opts = o; info.c = best.c;
end

function [best, nst] = bfgs_solve(A1t, A1p, A2t, A2p, S, g, keep, w1, w2, o, rs)
% Method 'bfgs': minimize J(c) = (w1 || |B1 c| - a1 ||^2 + w2 || |B2 c| - a2 ||^2)/2 over the
% kept coefficients with a quasi-Newton method (fminunc) and the exact gradient
% dJ/dconj(c) = sum_s w_s B_s' (B_s c - a_s .* exp(j angle(B_s c))) / 2.
% B1, B2 are small dense matrices (grid points x kept coefficients) built once.
% Starts: 1 = spectral initialization (leading eigenvector of sum_s w_s B_s' diag(a_s.^2) B_s),
% then o.Restarts - 1 random complex starts. Lowest J wins.
idx = find(keep); nc = numel(idx); nr = 2*numel(A1t);
B1 = zeros(nr, nc); B2 = B1;
for j = 1:nc
    e = zeros(numel(S.n), 1); e(idx(j)) = 1;
    [a, b] = swe_synth(S, e); B1(:, j) = [a(:); b(:)];
    B2(:, j) = B1(:, j) * g(idx(j));
end
a1 = [A1t(:); A1p(:)]; a2 = [A2t(:); A2p(:)];
M = w1*(B1' * (a1.^2 .* B1)) + w2*(B2' * (a2.^2 .* B2)); M = (M + M')/2;
[V, ~] = eigs(M, 1);
v = V(:, 1); sc = sqrt((w1*norm(a1)^2 + w2*norm(a2)^2) / (w1*norm(B1*v)^2 + w2*norm(B2*v)^2));
useTb = strcmp(o.Optimizer, 'fminunc') || (strcmp(o.Optimizer, 'auto') && has_optim());
if useTb
    opts = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'SpecifyObjectiveGradient', true, ...
        'Display', 'off', 'MaxIterations', o.MaxIter, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-14, ...
        'FunctionTolerance', 1e-14, 'MaxFunctionEvaluations', 5*o.MaxIter);
end
fun = @(x) cost(x, B1, B2, a1, a2, w1, w2);
best = struct('res', inf); nst = max(1, o.Restarts);
for t = 1:nst
    if t == 1 && ~isempty(o.InitCoef)
        c0 = o.InitCoef(idx);
    elseif t == 1 && ~isempty(o.InitPhase)
        c0 = swe_fit(S, A1t .* exp(1j*o.InitPhase{1}), A1p .* exp(1j*o.InitPhase{2})); c0 = c0(idx);
    elseif t == 1
        c0 = sc * v;
    else
        c0 = sc * (randn(rs, nc, 1) + 1j*randn(rs, nc, 1)) / sqrt(2) * norm(v) / sqrt(nc);
    end
    if useTb
        [x, J, ~, out] = fminunc(fun, [real(c0); imag(c0)], opts);
    else
        [x, J, out] = bfgs_min(fun, [real(c0); imag(c0)], o.MaxIter, 5*o.MaxIter, 1e-12, 1e-14, 1e-14);
    end
    if J < best.res
        c = zeros(numel(S.n), 1); c(idx) = x(1:nc) + 1j*x(nc+1:end);
        best = struct('res', J, 'c', c, 'iters', out.iterations, 'hist', J, 'restart', t);
    end
    if 10*log10(best.res) < o.StopRes_dB, nst = t; break; end
end
if o.Polish && ~has_optim()
    warning('twoscan_solve:polish', 'Polish needs the Optimization Toolbox (lsqnonlin); skipped');
elseif o.Polish
    x0 = [real(best.c(idx)); imag(best.c(idx))];
    lo = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'SpecifyObjectiveGradient', true, ...
        'Display', 'off', 'MaxIterations', 5*o.MaxIter, 'FunctionTolerance', 1e-20, 'StepTolerance', 1e-20, 'OptimalityTolerance', 1e-20);
    [x, rn, ~, ~, out] = lsqnonlin(@(x) lmres(x, B1, B2, a1, a2, w1, w2), x0, [], [], lo);
    if rn/2 < best.res
        c = zeros(numel(S.n), 1); c(idx) = x(1:nc) + 1j*x(nc+1:end);
        best.res = rn/2; best.c = c; best.iters = best.iters + out.iterations;
    end
end
end

function [J, gr] = cost(x, B1, B2, a1, a2, w1, w2)
nc = numel(x)/2; c = x(1:nc) + 1j*x(nc+1:end);
z1 = B1*c; z2 = B2*c;
d1 = abs(z1) - a1; d2 = abs(z2) - a2;
J = (w1*(d1'*d1) + w2*(d2'*d2)) / 2;
if nargout > 1
    u1 = z1 - a1 .* exp(1j*angle(z1)); u2 = z2 - a2 .* exp(1j*angle(z2));
    G = w1*(B1'*u1) + w2*(B2'*u2);            % = 2 dJ/dconj(c); real gradient [Re; Im]
    gr = [real(G); imag(G)];
end
end

function tf = has_optim()
tf = ~isempty(which('fminunc')) && license('test', 'Optimization_Toolbox');
end

function [x, f, out] = bfgs_min(fun, x, maxIter, maxEval, tolOpt, tolStep, tolFun)
%BFGS_MIN  Quasi-Newton (dense BFGS on the inverse Hessian) with a strong-Wolfe line search, a replacement for
%fminunc's 'quasi-newton' algorithm without the Optimization Toolbox (2026-10-19). fun returns [f, gradient].
%Stops when ||g||_inf <= tolOpt*max(1, ||g0||_inf), the step is below tolStep*(1 + ||x||), the decrease is below
%tolFun*(1 + |f|), or after maxIter iterations / maxEval evaluations. out.iterations, out.funcCount.
[f, g] = fun(x); nEval = 1; n = numel(x);
H = eye(n) / max(norm(g), eps);                       % first step: a short steepest-descent step
g0 = max(1, norm(g, inf)); it = 0;
while it < maxIter && nEval < maxEval
    it = it + 1;
    if norm(g, inf) <= tolOpt * g0, break; end
    d = -H * g;
    if g' * d >= 0, H = eye(n) / max(norm(g), eps); d = -H * g; end   % not a descent direction: reset
    [t, fn, gn, ne] = wolfe(fun, x, f, g, d);
    nEval = nEval + ne;
    if t == 0, break; end
    s = t * d; x = x + s; y = gn - g;
    df = f - fn; f = fn; g = gn;
    sy = s' * y;
    if sy > 1e-12 * norm(s) * norm(y)                  % BFGS update of the inverse Hessian
        if it == 1, H = (sy / (y' * y)) * eye(n); end   % initial scaling (Shanno-Phua)
        Hy = H * y; rho = 1 / sy;
        H = H + (rho^2 * (y' * Hy) + rho) * (s * s') - rho * (Hy * s' + s * Hy');
    end
    if norm(s) <= tolStep * (1 + norm(x)) || df <= tolFun * (1 + abs(f)), break; end
end
out.iterations = it; out.funcCount = nEval;
end

function [t, f, g, ne] = wolfe(fun, x, f0, g0, d)
% strong-Wolfe line search (c1 1e-4, c2 0.9) with bracketing and cubic/bisection zoom (Nocedal & Wright 3.5-3.6)
c1 = 1e-4; c2 = 0.9; dg0 = g0' * d; ne = 0;
tPrev = 0; fPrev = f0; dgPrev = dg0; t = 1;
for k = 1:25
    [f, g] = fun(x + t*d); ne = ne + 1; dg = g' * d;
    if f > f0 + c1*t*dg0 || (k > 1 && f >= fPrev)
        [t, f, g, n2] = zoom(fun, x, f0, dg0, d, tPrev, fPrev, dgPrev, t, f, dg, c1, c2); ne = ne + n2; return
    end
    if abs(dg) <= -c2*dg0, return, end
    if dg >= 0
        [t, f, g, n2] = zoom(fun, x, f0, dg0, d, t, f, dg, tPrev, fPrev, dgPrev, c1, c2); ne = ne + n2; return
    end
    tPrev = t; fPrev = f; dgPrev = dg; t = 2*t;
end
end

function [t, f, g, ne] = zoom(fun, x, f0, dg0, d, tlo, flo, dglo, thi, fhi, dghi, c1, c2)
ne = 0;
for k = 1:30
    t = cubic_min(tlo, flo, dglo, thi, fhi, dghi);
    [f, g] = fun(x + t*d); ne = ne + 1; dg = g' * d;
    if f > f0 + c1*t*dg0 || f >= flo
        thi = t; fhi = f; dghi = dg;
    else
        if abs(dg) <= -c2*dg0, return, end
        if dg*(thi - tlo) >= 0, thi = tlo; fhi = flo; dghi = dglo; end
        tlo = t; flo = f; dglo = dg;
    end
    if abs(thi - tlo) < 1e-14 * max(1, abs(tlo)), break; end
end
if flo < f0, t = tlo; [f, g] = fun(x + t*d); ne = ne + 1; else, t = 0; f = f0; g = []; end
if isempty(g), [f, g] = deal(f0, zeros(size(x))); end
end

function t = cubic_min(a, fa, ga, b, fb, gb)
% minimizer of the cubic through (a, fa, ga) and (b, fb, gb), kept inside the interval (else bisection)
d1 = ga + gb - 3*(fa - fb)/(a - b);
s = d1^2 - ga*gb;
lo = min(a, b); hi = max(a, b);
if s >= 0
    d2 = sign(b - a) * sqrt(s);
    t = b - (b - a) * (gb + d2 - d1) / (gb - ga + 2*d2);
    if ~(t > lo + 0.1*(hi - lo) && t < hi - 0.1*(hi - lo)), t = (a + b)/2; end
else
    t = (a + b)/2;
end
end

function [r, Jm] = lmres(x, B1, B2, a1, a2, w1, w2)
% residual vector sqrt(w_s) (|B_s c| - a_s) and its Jacobian w.r.t. [Re c; Im c]
nc = numel(x)/2; c = x(1:nc) + 1j*x(nc+1:end);
z1 = B1*c; z2 = B2*c;
r = [sqrt(w1)*(abs(z1) - a1); sqrt(w2)*(abs(z2) - a2)];
if nargout > 1
    u1 = conj(z1 ./ max(abs(z1), eps)); u2 = conj(z2 ./ max(abs(z2), eps));
    G1 = u1 .* B1; G2 = u2 .* B2;                 % d|z|/dRe c = Re(G), d|z|/dIm c = -Im(G)
    Jm = [sqrt(w1)*[real(G1), -imag(G1)]; sqrt(w2)*[real(G2), -imag(G2)]];
end
end
