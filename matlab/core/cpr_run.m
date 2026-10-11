function r = cpr_run(c, seed)
%CPR_RUN  One trial: pattern -> readings -> noise -> recovery -> score.
%   r = cpr_run(cpr_cfg('K', 8), 3)
% Everything random (probe/tone draws, cable lengths, noise, solver
% initializations) follows from rng(seed) in make_readings.m.
cpr_setup();                                       % puts GAMPmatlab (main, phase) on the path
[X, f] = pattern_crlh(c.nw, c.nd, c.f0idx);
x = X(:);
S = make_readings(c, f, seed);                     % calls rng(seed)
z = S.Phi * x;
wvar = 10^(-c.snrdB/10) * mean(abs(z).^2);
w = sqrt(wvar/2) * (randn(S.m, 2) * [1; 1i]);
PhiSolver = S.Phi; if isfield(S, 'PhiA'), PhiSolver = S.PhiA; end   % calibration error: solver assumes calibrated phases
A = SparseDftLinTrans(PhiSolver, c.nw, c.nh, c.nd, cpr_basis(c.basis, f, 0), c.sqMode);
t0 = tic;
switch c.meas
    case 'complex'
        y = z + w;
        o = struct();
        if c.gauge && c.nd > 1, o.gauge = @(e) cpr_gauge(e, A, 'delay', f, c.basis, S.rowTone); end
        [xh, info] = cpr_gamp_cs(y, A, o);
    case 'phaseless'
        y = abs(z + w);
        o = struct('minTry', c.minTry, 'maxTry', c.maxTry, 'nitEM', c.nitEM, 'SNRdB', c.solverSNRdB);
        if isfield(c, 'nitGAMP'), o.nitGAMP = c.nitGAMP; end
        if isfield(c, 'stopDB'), o.nresStopdB = c.stopDB; end
        o.reset = @() A.setM(cpr_basis(c.basis, f, 0));
        if c.gauge && c.nd > 1
            if c.K == 1, mode = 'perlambda'; else, mode = 'delay'; end
            o.gauge = @(e) cpr_gauge(e, A, mode, f, c.basis, S.rowTone);
        end
        if c.logRestarts, o.scoreFn = @(xp) cpr_score(X, xp, c.nd); end
        [xh, info] = feval(c.solver, y, A, o);
    otherwise
        error('cpr_run: unknown meas %s', c.meas);
end
r.cfg = c; r.seed = seed; r.time = toc(t0);
r.score = cpr_score(X, xh, c.nd);
r.info = info; r.R = S.R; r.m = S.m; r.n = S.n;
if isempty(xh), r.Xh = single([]); else, r.Xh = single(reshape(xh, size(X))); end
end
