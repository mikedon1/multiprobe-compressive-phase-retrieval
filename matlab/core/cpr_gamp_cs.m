function [xhat, info] = cpr_gamp_cs(y, A, o)
%CPR_GAMP_CS  Compressive sensing from COMPLEX measurements y = H*s + w with
%GAMP and the same EM-tuned Bernoulli-Gaussian prior as PR-GAMP (so complex
%and phaseless results differ only in the output channel). The noise
%variance is EM-tuned too. Linear GAMP needs no restarts.
%
%   o.gauge = @(estFin) ... : applied between GAMP rounds (delay re-choice,
%             cpr_gauge 'delay'); o.nRounds (default 10) rounds max.
% Returns the PATTERN estimate xhat = A.syn(shat) (A's basis may change).
if ~isfield(o, 'gauge'), o.gauge = []; end
if ~isfield(o, 'nRounds'), o.nRounds = 10; end
if ~isfield(o, 'damped'), o.damped = true; end
[m, n] = A.size();
k = 1; while (m > 2*k*log2(n/k)) && (k < n), k = k+1; end
rho = min(k/n, 0.9);
SNRinit = 1e3;                                          % 30 dB initial guess, then EM
wvar0 = mean(abs(y).^2) / (SNRinit + 1);
xvar0 = norm(y)^2 * SNRinit/(SNRinit+1) / rho / sum(A.multSq(ones(n,1)));
inner = CAwgnEstimIn(0, xvar0, false, 'mean0Tune', false, 'autoTune', true);
inputEst = SparseScaEstim(inner, rho, 0, 'autoTune', true);
outputEst = CAwgnEstimOut(y, wvar0, false, 'autoTune', true);
opt = GampOpt(); opt.legacyOut = 0; opt.nit = 300; opt.tol = 1e-7;
opt.adaptStep = 1; opt.stepWindow = 20; opt.uniformVariance = 0;
if o.damped   % structured (rotating) matrices: damp like PR-GAMP
    opt.step = 0.25; opt.stepMax = 0.5; opt.stepMin = 0.05; opt.adaptStepBethe = 1;
end
t0 = tic; xold = inf(n, 1); tau = 0;
for r = 1:o.nRounds
    [estFin, optFin] = gampEst(inputEst, outputEst, A, opt);
    xnew = A.syn(estFin.xhat);
    if ~isempty(o.gauge)
        estFin = o.gauge(estFin);
        if isfield(estFin, 'tau'), tau = estFin.tau; end
    end
    if norm(xnew - xold)/norm(xnew) < 1e-4, break; end
    xold = xnew;
    opt = optFin.warmStart(estFin, 'tol', 1e-7);
    if ~isempty(o.gauge), opt.valIn0 = -inf; opt.valOpt0 = []; end
end
xhat = A.syn(estFin.xhat);
info.rounds = r; info.time = toc(t0); info.tau = tau;
info.nres_dB = 20*log10(norm(A.Phi*xhat - y)/norm(y));
end
