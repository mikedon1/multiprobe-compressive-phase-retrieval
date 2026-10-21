function [Phi, d, info, g] = cpr_selfcal(ent, x, y, lam, d, nGN, sigL, opts)
%CPR_SELFCAL  Self-calibration step: per-port cable corrections from the current pattern.
%   [Phi, d, info] = cpr_selfcal(ent, x, y, lam, d, nGN, sigL)
%   [Phi, d, info, g] = cpr_selfcal(ent, x, y, lam, d, nGN, sigL, opts)
% Model: the solver's matrix entry e (row ent.i, column ent.j, port ent.port, tone ent.tone, assumed
% value ent.v) becomes ent.v * exp(g(port) + 1i*2*pi*d(port)/lam(tone)): each port's calibrated cable is off by a
% fixed length d (cm), so its phase error grows with frequency, and optionally by a fixed log-amplitude g. Given the
% pattern estimate x (in the span of the basis) and the amplitudes y, nGN damped Gauss-Newton steps minimize
% || |Phi(d, g) x| - y ||^2, starting from the current d and g. A length common to all ports is a delay that the
% pattern absorbs, and a gain common to all ports is the pattern's scale, so both means are held at zero.
% Optional sigL (cm): a zero-mean Gaussian prior of that rms on each length (MAP), weighted by the current
% data-residual variance; it keeps ports that appear in few readings from absorbing the pattern.
% opts (2026-10-19, all optional):
%   model  'length' (default): one length per port | 'drift': one parameter a, d = a*(Lcal - mean(Lcal)) (thermal
%          drift: every cable longer by the same fraction; needs opts.Lcal, the calibrated lengths per port, cm)
%   gain   true: also one log-amplitude g per port (prior opts.sigG, natural log; NaN/absent = none)
%   g      current log-amplitudes (default zeros)
%   projectAfter  true (default): the step used for the paper's results (solve the damped normal equations, then
%          subtract the mean of the step); false: the constrained Gauss-Newton step, solved on an orthonormal basis Q
%          of zero-mean corrections, step = -Q (Q' H Q + mu I) \ Q' grad (round-8 review). At sigma 5 deg, seed 1, the
%          two differ by at most 0.15 dB (IC -30.31/-30.46, two-pass IC -33.39/-33.45, IS -32.31/-32.28; diag_selfcal)
% Returns the re-phased matrix, the new d (and g), and info (residual before/after, dB).
if nargin < 6 || isempty(nGN), nGN = 5; end
if nargin < 7, sigL = NaN; end
if nargin < 8, opts = struct(); end
model = 'length'; if isfield(opts, 'model'), model = opts.model; end
useG = isfield(opts, 'gain') && opts.gain;
projectAfter = ~isfield(opts, 'projectAfter') || opts.projectAfter;   % paper's step (solve, then subtract the mean)
lam = lam(:); nP = ent.nPorts; m = ent.m;
Q = null(ones(1, nP));                               % orthonormal basis of zero-mean per-port corrections
g = zeros(nP, 1); if isfield(opts, 'g') && ~isempty(opts.g), g = opts.g(:); end
k = 1i*2*pi ./ lam(ent.tone);                        % d(phase)/d(length) per entry
k = k(:);
xe = x(ent.j);
if strcmp(model, 'drift')
    Lc = opts.Lcal(:) - mean(opts.Lcal(:));
    a0 = (Lc' * d) / max(Lc' * Lc, eps); d = a0 * Lc;   % project the current d on the drift direction
end
fwd = @(dd, gg) ent.v .* exp(gg(ent.port) + k .* dd(ent.port)) .* xe;
data = @(dd, gg) norm(abs(accumarray(ent.i, fwd(dd, gg), [m 1])) - y)^2;
betaL = 0; betaG = 0; nv = data(d, g)/m;              % noise variance estimate from the current residual
if isfinite(sigL) && sigL > 0 && ~strcmp(model, 'drift'), betaL = nv / sigL^2; end
if useG && isfield(opts, 'sigG') && isfinite(opts.sigG) && opts.sigG > 0, betaG = nv / opts.sigG^2; end
cost = @(dd, gg) data(dd, gg) + betaL*norm(dd)^2 + betaG*norm(gg)^2;
c0 = cost(d, g); c = c0;
for it = 1:nGN
    a = fwd(d, g);
    z = accumarray(ent.i, a, [m 1]);
    r = abs(z) - y;
    u = conj(z) ./ max(abs(z), eps);
    Jd = sparse(ent.i, ent.port, real(u(ent.i) .* k .* a), m, nP);   % d|z_i| / d d_p
    switch model
        case 'drift', J = Jd * Lc; H = J'*J; gr = J'*r;            % one parameter a
        otherwise,    J = Jd;      H = J'*J + betaL*speye(nP); gr = J'*r + betaL*d;
    end
    if useG
        Jg = sparse(ent.i, ent.port, real(u(ent.i) .* a), m, nP); % d|z_i| / d g_p
        H = [H, J'*Jg; Jg'*J, Jg'*Jg + betaG*speye(nP)]; gr = [gr; Jg'*r + betaG*g]; %#ok<AGROW>
    end
    nq = size(H, 1);
    mu = 1e-6 * full(sum(diag(H))) / nq + eps;
    if strcmp(model, 'drift')                       % one parameter; Lc already has zero mean
        step = -((H + mu*speye(nq)) \ gr); sd = step(1) * Lc; ns = 1; sg = zeros(nP, 1);
        if useG, sg = step(ns+1:end); sg = sg - mean(sg); end
    elseif projectAfter                             % paper's step: solve, then remove the mean
        step = -((H + mu*speye(nq)) \ gr); sd = step(1:nP) - mean(step(1:nP)); ns = nP;
        sg = zeros(nP, 1); if useG, sg = step(ns+1:end); sg = sg - mean(sg); end
    else                                            % constrained Gauss-Newton step in the zero-mean subspace (2026-10-21)
        Qb = Q; if useG, Qb = blkdiag(Q, Q); end
        Hq = full(Qb' * H * Qb); Hq = (Hq + Hq')/2;
        stepq = -((Hq + mu*eye(size(Hq, 1))) \ (Qb' * gr)); step = Qb * stepq;
        sd = step(1:nP); sg = zeros(nP, 1); if useG, sg = step(nP+1:end); end
    end
    t = 1; improved = false;
    while t > 1/64
        dn = d + t*sd; gn = g + t*sg; cn = cost(dn, gn);
        if cn < c, d = dn; g = gn; c = cn; improved = true; break; end
        t = t/2;
    end
    if ~improved || (c0 - c) < 1e-9*c0, break; end
end
Phi = sparse(ent.i, ent.j, ent.v .* exp(g(ent.port) + k .* d(ent.port)), m, ent.n);
info.res0_dB = 10*log10(c0/norm(y)^2); info.res_dB = 10*log10(data(d, g)/norm(y)^2);
end
