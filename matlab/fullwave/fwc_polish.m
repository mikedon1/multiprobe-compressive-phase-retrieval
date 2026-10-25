function [ah, info] = fwc_polish(A, y, ah0, varargin)
%FWC_POLISH  Support-growing magnitude least-squares refinement of a PR-GAMP estimate (Section VII only).
%   [ah, info] = fwc_polish(A, y, ah0)
% A    CalVshLinTrans operator;  y  measured amplitudes;  ah0  PR-GAMP pattern estimate (A.syn order, p*nf x 1 or p x nf)
% PR-GAMP accepts a restart at a residual of -(SNR_s + 2) dB, well above the noise level when the pattern is only weakly
% sparse (2026-10-23 study). This step minimizes 0.5 || |A s| - y ||^2 over the coefficients s of a support that starts
% at the PR-GAMP support (|s| > 1e-3 max) and grows by 'Step' coefficients per round (the ones with the largest
% gradient), each round solved with L-BFGS. No truth and no noise level are used. The round is chosen by a held-out set:
% a random 'HoldOut' fraction of the readings is not fitted, and the round with the smallest residual on it is kept.
% Options: 'Step' 250, 'Rounds' 14, 'HoldOut' 0.1, 'Iters' 200 (L-BFGS iterations per round), 'Seed' 1, 'Refit' true (refit the chosen support on all readings; false: keep the fit on the 90%).
% info: round (chosen), heldoutDB / trainDB (per round, dB), support (size per round).
p = inputParser;
p.addParameter('Step', 250); p.addParameter('Rounds', 14); p.addParameter('HoldOut', 0.1);
p.addParameter('Iters', 200); p.addParameter('Seed', 1); p.addParameter('Refit', true);
p.parse(varargin{:}); o = p.Results;
m = numel(y); y = y(:);
rs = RandStream('mt19937ar', 'Seed', 500 + o.Seed);
wt = double(rand(rs, m, 1) > o.HoldOut); ho = wt == 0;
s0 = A.ana(ah0(:)); mask = abs(s0) > 1e-3*max(abs(s0)); s = s0 .* mask;
hr = zeros(1, o.Rounds); tr = hr; sup = hr; best = inf; sBest = s; ib = 1; mBest = mask;
for it = 1:o.Rounds
    s = lbfgs_mag(A, y, s, mask, o.Iters, wt);
    r = abs(A.mult(s)) - y;
    tr(it) = 10*log10(sum((wt.*r).^2) / sum((wt.*y).^2)); hr(it) = 10*log10(sum(r(ho).^2) / sum(y(ho).^2));
    sup(it) = nnz(mask);
    if hr(it) < best, best = hr(it); sBest = s; ib = it; mBest = mask; end
    [~, g] = fg(A, y, s, true(size(s)), wt); g(mask) = 0; [~, o2] = sort(abs(g), 'descend'); mask(o2(1:o.Step)) = true;
end
if o.Refit, sBest = lbfgs_mag(A, y, sBest, mBest, o.Iters, ones(m, 1)); end   % final fit of the chosen support on ALL readings
ah = A.syn(sBest);
info = struct('round', ib, 'heldoutDB', hr, 'trainDB', tr, 'support', sup);
end

function s = lbfgs_mag(A, y, s, mask, nit, wt)
mem = 10; S = {}; Y = {};
[f, g] = fg(A, y, s, mask, wt);
for it = 1:nit
    q = g; al = zeros(numel(S), 1);
    for k = numel(S):-1:1, rho = 1/real(Y{k}'*S{k}); al(k) = rho*real(S{k}'*q); q = q - al(k)*Y{k}; end
    if ~isempty(S), q = q * real(S{end}'*Y{end})/real(Y{end}'*Y{end}); end
    for k = 1:numel(S), rho = 1/real(Y{k}'*S{k}); be = rho*real(Y{k}'*q); q = q + S{k}*(al(k)-be); end
    d = -q; if real(d'*g) >= 0, d = -g; S = {}; Y = {}; end
    t = 1; ok = false;
    for ls = 1:30
        sn = s + t*d; [fn, gn] = fg(A, y, sn, mask, wt);
        if fn <= f + 1e-4*t*real(g'*d), ok = true; break; end
        t = t/2;
    end
    if ~ok, break; end
    S{end+1} = sn - s; Y{end+1} = gn - g; if numel(S) > mem, S(1) = []; Y(1) = []; end %#ok<AGROW>
    s = sn; f = fn; g = gn;
end
end

function [f, g] = fg(A, y, s, mask, wt)
z = A.mult(s); az = abs(z); r = wt .* (az - y); f = 0.5*(r'*r);
g = A.multTr(r .* z ./ max(az, eps)); g = g .* mask;
end
