function S = swe_setup(theta, phi, N)
%SWE_SETUP  Vector spherical-wave expansion on a theta x phi tensor grid (two-sphere IFT).
%
%   S = swe_setup(theta, phi, N)   theta: nt x 1 in [0, pi]; phi: np x 1 equispaced over a
%                                  full circle (np >= 2N+1); N: highest mode degree.
%
% Tangential field on a sphere of radius r:
%   E_t = sum_{n,m} cTE_nm(r) X_nm  +  cTM_nm(r) (rhat x X_nm),      n = 1..N, |m| <= n
% with X_nm = L Y_nm / sqrt(n(n+1)) (vector spherical harmonic):
%   X      : theta comp  -m Ybar/sin(theta),   phi comp  -1i dYbar/dtheta
%   rhat x X: theta comp  1i dYbar/dtheta,     phi comp  -m Ybar/sin(theta)
% (all times exp(1i m phi)/sqrt(n(n+1)); Ybar = fully normalized Legendre P_n^|m|(cos theta)
% / sqrt(2 pi)). The sign convention for negative m is irrelevant here: the same basis is used
% to fit and to synthesize, and propagation acts per n (swe_radial.m).
% Radial dependence (exp(+j w t), outgoing): cTE ~ h_n^(2)(kr), cTM ~ (1/kr) d/dkr[kr h_n^(2)].
%
% S.A{mi}   : (2 nt) x nc basis for m = S.m(mi) (rows: theta comps then phi comps)
% S.Pi{mi}  : sin(theta)-weighted least-squares pseudo-inverse of S.A{mi}
% S.n, S.type, S.mi : degree, 1 (TE) / 2 (TM), and m-block of every coefficient (vector form)
% S.Ephi    : np x (2N+1) matrix exp(1i phi m)
theta = theta(:); phi = phi(:); nt = numel(theta); np = numel(phi);
if np < 2*N + 1, error('swe_setup: need np >= 2N+1 phi samples'); end
th = min(max(theta, 1e-4), pi - 1e-4);                     % keep the pole rows finite (and theta-h > 0)
h = 1e-5;                                                  % theta step for d/dtheta
w = sqrt(sin(th)) + 1e-3;                                  % LS weight ~ sqrt(dOmega)
S.m = -N:N; S.N = N; S.theta = theta; S.phi = phi;
S.Ephi = exp(1i * phi * S.m);
S.A = cell(1, 2*N + 1); S.Pi = S.A; S.n = []; S.type = []; S.mi = [];
P0 = cell(1, N); Pp = P0; Pm = P0;
for n = 1:N
    P0{n} = legendre(n, cos(th), 'norm') / sqrt(2*pi);     % (n+1) x nt, rows m = 0..n
    Pp{n} = legendre(n, cos(th + h), 'norm') / sqrt(2*pi);
    Pm{n} = legendre(n, cos(th - h), 'norm') / sqrt(2*pi);
end
for mi = 1:numel(S.m)
    m = S.m(mi); ns = max(1, abs(m)):N; nc = numel(ns);
    A = zeros(2*nt, 2*nc);
    for j = 1:nc
        n = ns(j); c = 1/sqrt(n*(n + 1));
        Y = P0{n}(abs(m) + 1, :).';
        dY = (Pp{n}(abs(m) + 1, :).' - Pm{n}(abs(m) + 1, :).') / (2*h);
        A(:, j)      = c * [-m*Y./sin(th); -1i*dY];        % TE: X
        A(:, nc + j) = c * [1i*dY; -m*Y./sin(th)];          % TM: rhat x X
    end
    S.A{mi} = A;
    W = [w; w];
    S.Pi{mi} = pinv(W .* A) .* W.';
    S.n = [S.n, ns, ns]; S.type = [S.type, ones(1, nc), 2*ones(1, nc)];
    S.mi = [S.mi, mi*ones(1, 2*nc)];
end
S.n = S.n(:); S.type = S.type(:); S.mi = S.mi(:);
end
