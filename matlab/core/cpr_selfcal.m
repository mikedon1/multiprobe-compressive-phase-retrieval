function [Phi, d, info] = cpr_selfcal(ent, x, y, lam, d, nGN, sigL)
%CPR_SELFCAL  Self-calibration step: per-port cable length corrections from the current pattern.
%   [Phi, d, info] = cpr_selfcal(ent, x, y, lam, d, nGN)
% Model: the solver's matrix entry e (row ent.i, column ent.j, port ent.port, tone ent.tone, assumed
% value ent.v) becomes ent.v * exp(1i*2*pi*d(port)/lam(tone)): each port's calibrated cable is off by a
% fixed length d (cm), so its phase error grows with frequency. Given the pattern estimate x (in the span
% of the basis) and the amplitudes y, nGN damped Gauss-Newton steps on d minimize || |Phi(d) x| - y ||^2,
% starting from the current d. A length common to all ports is a delay that the pattern absorbs, so the
% mean of d is held at zero. Optional sigL (cm): a zero-mean Gaussian prior of that rms on each length
% (MAP), weighted by the current data-residual variance; it keeps ports that appear in few readings from
% absorbing the pattern. Returns the re-phased matrix, the new d, and info (residual before/after, dB).
if nargin < 6 || isempty(nGN), nGN = 5; end
if nargin < 7, sigL = NaN; end
lam = lam(:); nP = ent.nPorts; m = ent.m;
k = 1i*2*pi ./ lam(ent.tone);                        % d(phase)/d(length) per entry
xe = x(ent.j);
data = @(dd) norm(abs(accumarray(ent.i, ent.v .* exp(k .* dd(ent.port)) .* xe, [m 1])) - y)^2;
beta = 0;
if isfinite(sigL) && sigL > 0, beta = (data(d)/m) / sigL^2; end   % noise variance / prior variance
cost = @(dd) data(dd) + beta*norm(dd)^2;
c0 = cost(d); c = c0;
for it = 1:nGN
    a = ent.v .* exp(k .* d(ent.port)) .* xe;
    z = accumarray(ent.i, a, [m 1]);
    r = abs(z) - y;
    u = conj(z) ./ max(abs(z), eps);
    Jm = sparse(ent.i, ent.port, real(u(ent.i) .* k .* a), m, nP);   % d|z_i| / d d_p
    H = Jm' * Jm + beta*speye(nP); g = Jm' * r + beta*d;
    mu = 1e-6 * full(sum(diag(H))) / nP + eps;
    step = -((H + mu*speye(nP)) \ g);
    step = step - mean(step);                         % gauge: no common length change
    t = 1; improved = false;
    while t > 1/64
        dn = d + t*step; cn = cost(dn);
        if cn < c, d = dn; c = cn; improved = true; break; end
        t = t/2;
    end
    if ~improved || (c0 - c) < 1e-9*c0, break; end
end
Phi = sparse(ent.i, ent.j, ent.v .* exp(k .* d(ent.port)), m, ent.n);
info.res0_dB = 10*log10(c0/norm(y)^2); info.res_dB = 10*log10(data(d)/norm(y)^2);
end
