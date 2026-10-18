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
%
% info: numTry (restarts run), nres_dB (best residual), time (s), log (per restart).

def = struct('SNRdB', 10, 'stopDB', NaN, 'minTry', 1, 'maxTry', 10, 'nitEM', 20, ...
    'nitGAMP', 200, 'gauge', [], 'reset', [], 'scoreFn', []);
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
P = init_params(yAmp, A, m, n, snrS);

% ----- prior and GAMP settings shared by all restarts
coefEst = CAwgnEstimIn(zeros(n,1), P.v0, false, 'mean0Tune', false, 'autoTune', true);
priorEst = SparseScaEstim(coefEst, P.rho0, 0, 'autoTune', true);
gopt = gamp_settings(opts.nitGAMP);
gain = 100;                                   % initial-estimate gain, raised after a failure

best = struct('nres', inf, 'x', []);
info.log = struct('nres_dB', {}, 'time', {}, 'score', {});
tStart = tic; tryCount = 0;
for t = 1:nTry
    tryCount = t; tRestart = tic;
    try
        if ~isempty(opts.reset), opts.reset(); end
        [s, nres, failed] = one_restart(yAmp, A, P, coefEst, priorEst, gopt, gain, opts);
        if failed
            gain = 2*gain;                    % numerical failure: larger start, shorter step window
            gopt.stepWindow = ceil(0.75*gopt.stepWindow);
        end
        xPat = A.syn(s);
        entry.nres_dB = nres; entry.time = toc(tRestart); entry.score = [];
        if ~isempty(opts.scoreFn), entry.score = opts.scoreFn(xPat); end
        info.log(t) = entry;
        if nres < best.nres, best.nres = nres; best.x = xPat; end
        if t >= opts.minTry && best.nres < target, break; end
    catch err
        warning('cpr_prgamp:restart', 'restart %d failed: %s', t, err.message);
    end
end
xhat = best.x;
info.numTry = tryCount; info.nres_dB = best.nres; info.time = toc(tStart);
end

% -------------------------------------------------------------------------------------------
function P = init_params(yAmp, A, m, n, snrS)
% Sparsity rate: smallest k with m <= 2 k log2(n/k), capped at 0.9 n. Coefficient variance:
% measured energy, scaled by the SNR, shared over rho*n active coefficients and ||A||_F^2.
k = 1;
while k < n && m > 2*k*log2(n/k), k = k + 1; end
P.rho0 = min(k/n, 0.9);
snrLin = 10^(snrS/10);
P.v0 = (norm(yAmp)^2*snrLin/(snrLin+1)/P.rho0/sum(A.multSq(ones(n,1))))*ones(n,1);
P.wvar0 = mean(abs(yAmp).^2)/(10 + 1);        % noise variance starts at 10 dB SNR
P.tolMax = 1e-4;
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

function [s, nres, failed] = one_restart(yAmp, A, P, coefEst, priorEst, gopt, gain, opts)
% One restart: random start, then alternate GAMP runs and EM updates.
n = numel(P.v0); m = numel(yAmp);
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
    nres = 20*log10(norm(abs(A.mult(est.xhat)) - yAmp)/norm(yAmp));
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
