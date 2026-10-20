function [xhat, info] = cpr_prgamp(yAmp, A, opts)
%CPR_PRGAMP  Phase retrieval with GAMP: restarts, EM tuning, and early stopping.
%
%   [xhat, info] = cpr_prgamp(yAmp, A, opts)
%
% Recovers sparse coefficients s from amplitudes yAmp = |A*s + w| and returns the PATTERN
% estimate xhat = A.syn(s_hat) of the restart with the lowest amplitude residual. Written for
% this project (BSD 3-Clause Clear license) from the procedure stated in the paper (Section IV and
% Appendix B). It calls GAMPmatlab (gampEst, GampOpt, CAwgnEstimIn, SparseScaEstim,
% ncCAwgnEstimOut), which must be on the MATLAB path; GAMPmatlab itself is not part of this
% code.
%
% Procedure
%   Prior: Bernoulli-Gaussian on s (sparsity rate rho, variance v), learned by EM.
%   Noise: complex Gaussian before detection, variance learned by EM.
%   Each restart starts GAMP from a random complex Gaussian estimate and alternates GAMP runs
%   with EM updates. A restart ends when the estimate stops changing or after nEM updates.
%   The solver stops at the first restart whose residual R = 20 log10(|| |A s| - y || / ||y||)
%   is below the target, or keeps the best of maxTry restarts.
%
% opts (all optional)
%   SNRdB      solver setting SNR_s; residual target -(SNR_s + 2) dB          (10)
%   stopDB     explicit residual target in dB, overrides the SNR_s rule     (NaN = unused)
%   minTry     restarts always run before early stopping is allowed          (1)
%   maxTry     maximum number of restarts                                    (10)
%   nitEM      EM updates per restart (nitEM + 1 GAMP runs)                  (20)
%   nitGAMP    GAMP iterations per run                                       (200)
%   gauge      @(est) re-expression of the estimate after each EM update     ([])
%   reset      @() called at the start of every restart                     ([])
%   scoreFn    @(xpattern) score logged for every restart                    ([])
%   verbose    print one progress line per restart                           (false)
%   groupPrior structured sparsity rate (VSH basis, A.npol = 2, test 2026-10-18):  ('none')
%              'pair'    one rate per (l, m, frequency-basis index): TE and TM share it
%              'pairsep' rate = (l, m) group over all frequency indices x frequency-index factor
%              'modesep' rate = mode group over all frequency indices x frequency-index factor
%   sparsePrior false: band-limit-only ablation, Gaussian prior on every coefficient          (true)
%
% info: numTry (restarts run), nres_dB (best residual), time (s), log (per restart).

def = struct('SNRdB', 10, 'stopDB', NaN, 'minTry', 1, 'maxTry', 10, 'nitEM', 20, ...
    'nitGAMP', 200, 'gauge', [], 'reset', [], 'scoreFn', [], 'verbose', false, 'groupPrior', 'none', ...
    'sparsePrior', true);
if nargin < 3, opts = struct(); end
for f = fieldnames(def).'
    if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
end
if isfield(opts, 'nresStopdB'), opts.stopDB = opts.nresStopdB; end    % older option name

[m, n] = A.size();
nTry = max(opts.minTry, opts.maxTry);
snrS = min(100, opts.SNRdB);
target = -(snrS + 2);
if ~isnan(opts.stopDB), target = opts.stopDB; end

% ----- initial prior and noise parameters (Appendix B)
P = init_params(yAmp, A, m, n, snrS, opts.sparsePrior);

% ----- prior and GAMP settings shared by all restarts
coefEst = CAwgnEstimIn(zeros(n,1), P.v0, false, 'mean0Tune', false, 'autoTune', true);
if ~opts.sparsePrior
    % band-limit-only ablation: rate 1 makes the Bernoulli-Gaussian prior a plain Gaussian on every
    % coefficient (log(1 - p1) = -inf gives activity 1, and the EM update keeps p1 = 1); variance learned
    priorEst = SparseScaEstim(coefEst, 1, 0, 'autoTune', true);
elseif strcmp(opts.groupPrior, 'none')
    priorEst = SparseScaEstim(coefEst, P.rho0, 0, 'autoTune', true);
else
    [grp, col, gmode] = group_labels(A, opts.groupPrior);
    priorEst = GroupSparseScaEstim(coefEst, P.rho0, 0, grp, gmode, 'autoTune', true);
    priorEst.cols = col;
end
gopt = gamp_settings(opts.nitGAMP);
gain = 100;                                   % initial-estimate gain, raised after a failure

best = struct('nres', inf, 'x', []);
info.log = struct('nres_dB', {}, 'time', {}, 'score', {}, 'nRuns', {}, 'nIts', {});
tStart = tic; tryCount = 0;
for t = 1:nTry
    tryCount = t; tRestart = tic;
    try
        if ~isempty(opts.reset), opts.reset(); end
        [s, nres, failed, nRuns, nIts] = one_restart(yAmp, A, P, coefEst, priorEst, gopt, gain, opts);
        if failed
            gain = 2*gain;                    % numerical failure: larger start, shorter step window
            gopt.stepWindow = ceil(0.75*gopt.stepWindow);
        end
        xPat = A.syn(s);
        entry.nres_dB = nres; entry.time = toc(tRestart); entry.score = [];
        if ~isempty(opts.scoreFn), entry.score = opts.scoreFn(xPat); end
        entry.nRuns = nRuns; entry.nIts = nIts;   % GAMP runs (EM updates + 1) and GAMP iterations of this restart (round 7)
        info.log(t) = entry;
        if nres < best.nres, best.nres = nres; best.x = xPat; end
        if opts.verbose
            fprintf('  restart %d/%d: residual %.1f dB (target %.1f), %.0f s  [%s]\n', ...
                t, nTry, nres, target, entry.time, datestr(now));
        end
        if t >= opts.minTry && best.nres < target, break; end
    catch err
        warning('cpr_prgamp:restart', 'restart %d failed: %s', t, err.message);
    end
end
xhat = best.x;
info.numTry = tryCount; info.nres_dB = best.nres; info.time = toc(tStart);
info.nRunsTotal = sum([info.log.nRuns]); info.nItsTotal = sum([info.log.nIts]);   % over the logged restarts
end

% -------------------------------------------------------------------------------------------
function P = init_params(yAmp, A, m, n, snrS, sparse)
% Sparsity rate: smallest k with m <= 2 k log2(n/k), capped at 0.9 n (1 without the sparsity prior).
% Coefficient variance: measured energy, scaled by the SNR, shared over rho*n active coefficients and ||A||_F^2.
k = 1;
while k < n && m > 2*k*log2(n/k), k = k + 1; end
P.rho0 = min(k/n, 0.9);
if ~sparse, P.rho0 = 1; end
snrLin = 10^(snrS/10);
P.v0 = (norm(yAmp)^2*snrLin/(snrLin+1)/P.rho0/sum(A.multSq(ones(n,1))))*ones(n,1);
P.wvar0 = mean(abs(yAmp).^2)/(10 + 1);        % noise variance starts at 10 dB SNR
P.tolMax = 1e-4;
end

function [grp, col, gmode] = group_labels(A, kind)
% Coefficient s = reshape(p, nf): column j of the VSH basis (vsh_basis, degree order; with npol = 2
% no column is dropped, so column j is mode j of vsh_eval up to the Gram-Schmidt mixing) times
% frequency-basis index f.
if ~(isa(A, 'VshLinTrans') || isa(A, 'CalVshLinTrans'))
    error('cpr_prgamp: groupPrior needs a VSH operator (VshLinTrans or CalVshLinTrans)');
end
[J, F] = ndgrid(1:A.p, 1:A.nf);
col = F(:);
if strcmp(kind, 'modesep')                     % groups = basis columns: any npol
    grp = J(:); gmode = 'sep'; return
end
if A.npol ~= 2 || A.p ~= 2*A.L*(A.L+2)          % (l, m) pairing needs the column = mode identity
    error('cpr_prgamp: groupPrior %s needs the dual-polarized VSH operator with all 2L(L+2) columns', kind);
end
[~, ~, modes] = vsh_eval(A.L, 1, 0);
[~, ~, lm] = unique(modes(:, 1:2), 'rows', 'stable');      % (l, m) index per column
nlm = max(lm);
switch kind
    case 'pair',    grp = lm(J(:)) + nlm*(F(:) - 1); gmode = 'avg';
    case 'pairsep', grp = lm(J(:));                  gmode = 'sep';
    case 'modesep', grp = J(:);                      gmode = 'sep';
    otherwise, error('cpr_prgamp: unknown groupPrior %s', kind);
end
end

function g = gamp_settings(nit)
% Adaptive damping with step at most 0.25, initial step window 50 (Appendix B).
g = GampOpt();
g.legacyOut = 0; g.nit = nit; g.histIntvl = 1;
g.adaptStep = 1; g.adaptStepBethe = 1;
g.step = 0.25; g.stepMax = 0.25; g.stepMin = 0.05;
g.stepWindow = 50; g.stepIncr = 1.1; g.stepDecr = 0.2;
g.varNorm = 1; g.pvarStep = 1; g.rvarStep = 0; g.uniformVariance = 0;
end

function [s, nres, failed, nRuns, nIts] = one_restart(yAmp, A, P, coefEst, priorEst, gopt, gain, opts)
% One restart: random start, then alternate GAMP runs and EM updates. nRuns, nIts: GAMP runs and iterations used.
n = numel(P.v0); m = numel(yAmp); nRuns = 0; nIts = 0;
start = sqrt(P.rho0*P.v0/2).*(randn(n,2)*[1;1i]);
gopt.xhat0 = gain*start;
gopt.xvar0 = mean(abs(gopt.xhat0).^2)*ones(n,1);
priorEst.p1 = P.rho0; coefEst.mean0 = zeros(n,1); coefEst.var0 = P.v0;   % fresh prior
wvar = P.wvar0; xPrev = inf(n,1); reexpressed = false; failed = false;
for it = 1:opts.nitEM + 1
    snrHat = 10*log10(mean(abs(yAmp).^2)/wvar);
    tolG = min(10^(-snrHat/10), P.tolMax);
    if it == 1                                % first run: prior and noise fixed
        run = gopt; run.tol = tolG; run.xvar0auto = 0;
        noiseEst = ncCAwgnEstimOut(yAmp, wvar*ones(m,1), 0, 0);
        priorEst.disableTune = true;
    else                                      % later runs: warm start, EM tuning on
        run = lastOpt.warmStart(est, 'tol', tolG); run.xvar0auto = false;
        if reexpressed, run.valIn0 = -inf; run.valOpt0 = []; end
        if it == 2, noiseEst = ncCAwgnEstimOut(yAmp, wvar*ones(m,1), 0, 1); end
        priorEst.disableTune = false;
    end
    [est, lastOpt] = gampEst(priorEst, noiseEst, A, run);
    nRuns = it; if isfield(est, 'nit'), nIts = nIts + est.nit; end
    nres =20*log10(norm(abs(A.mult(est.xhat)) - yAmp)/norm(yAmp));
    collapsed = norm(est.rhat)^2 <= 1e-6*P.rho0*sum(P.v0);
    if ~isfinite(nres) || collapsed
        failed = true; break;
    end
    reexpressed = false;
    if ~isempty(opts.gauge) && it <= opts.nitEM
        est = opts.gauge(est); reexpressed = true;
    end
    if it == opts.nitEM + 1 || norm(xPrev - est.xhat)/norm(est.xhat) < 10*tolG, break; end
    xPrev = est.xhat;
    wvar = mean(noiseEst.var0);
end
s = est.xhat;
end
