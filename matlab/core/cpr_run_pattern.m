function r = cpr_run_pattern(c, seed, X, f)
%CPR_RUN_PATTERN  cpr_run.m for an arbitrary pattern X (nw x nh x nd*npol) and frequencies f.
%Used by the two-scan comparison (near-field data). c.nw and c.nh are set from X, so the grid need
%not be nh = nw/2. With c.npol = 2 (dual polarization, Section VII), X stacks the polarizations
%along the third dimension ([pol 1: f_1..f_nd, pol 2: f_1..f_nd]); the basis is the same per
%polarization with one shared delay (cpr_basis), and readings combine ports of both
%(make_readings). Everything else (readings, noise, solver, gauge, scoring) is cpr_run.m.
% r.score is over all nd*npol slices; r.scoreVec is the vector score: one phase per wavelength
% shared by both polarizations ('perlambda' column) and one global phase ('global').
cpr_setup();                                       % puts GAMPmatlab (main, phase) on the path
if ~isfield(c, 'npol'), c.npol = 1; end
[c.nw, c.nh, ns] = size(X); c.nd = ns / c.npol; c.K = min(c.K, c.nd);
if c.nd ~= numel(f), error('cpr_run_pattern: size(X,3) must be npol*numel(f)'); end
x = X(:);
S = make_readings(c, f, seed);
z = S.Phi * x;
wvar = cpr_noisevar(c, z, x);
w = sqrt(wvar/2) * (randn(S.m, 2) * [1; 1i]);
PhiSolver = S.Phi; if isfield(S, 'PhiA'), PhiSolver = S.PhiA; end   % calibration error: solver assumes calibrated phases
A = cpr_operator(c, PhiSolver, c.nw, c.nh, ns, cpr_basis(c.basis, f, 0, c.npol));
t0 = tic;
switch c.meas
    case 'complex'
        o = struct();
        if c.gauge && c.nd > 1, o.gauge = @(e) cpr_gauge(e, A, 'delay', f, c.basis, S.rowTone); end
        [xh, info] = cpr_gamp_cs(z + w, A, o);
    case 'phaseless'
        o = struct('minTry', c.minTry, 'maxTry', c.maxTry, 'nitEM', c.nitEM, 'SNRdB', min(c.solverSNRdB, c.snrdB));
        if isfield(c, 'nitGAMP'), o.nitGAMP = c.nitGAMP; end
        if isfield(c, 'stopDB'), o.nresStopdB = c.stopDB; end
        if isfield(c, 'verbose'), o.verbose = c.verbose; end
        if isfield(c, 'groupPrior'), o.groupPrior = c.groupPrior; end
        if isfield(c, 'sparsePrior'), o.sparsePrior = c.sparsePrior; end
        o.reset = @() A.setM(cpr_basis(c.basis, f, 0, c.npol));
        if c.gauge && c.nd > 1
            if c.K == 1, mode = 'perlambda'; else, mode = 'delay'; end
            o.gauge = @(e) cpr_gauge(e, A, mode, f, c.basis, S.rowTone);
        end
        if isfield(c, 'logRestarts') && c.logRestarts, o.scoreFn = @(xp) vec_score(X, xp, c.nd, c.npol); end
        y = abs(z + w);
        selfCal = isfield(c, 'selfCal') && c.selfCal;
        if selfCal
            % self-calibration: after each EM update from c.selfCalStart on, re-estimate one cable length
            % correction per port from the current (in-span) pattern, re-phase the solver's matrix, then
            % apply the usual re-expression. Every restart starts from the assumed calibration.
            calD = zeros(S.ent.nPorts, 1); calG = calD; calCalls = 0; calLog = [];
            calOpts = struct('model', 'length', 'gain', false, 'Lcal', S.L(:));
            if isfield(c, 'selfCalModel'), calOpts.model = c.selfCalModel; end
            if isfield(c, 'selfCalGain'), calOpts.gain = c.selfCalGain; end
            if isfield(c, 'selfCalPriorGainDb'), calOpts.sigG = c.selfCalPriorGainDb/20*log(10); end
            gauge0 = []; if isfield(o, 'gauge'), gauge0 = o.gauge; end
            reset0 = o.reset;
            o.reset = @calReset;
            o.gauge = @calHook;
        end
        [xh, info] = feval(c.solver, y, A, o);
        if selfCal
            info.selfCal = cal_score(calD, S, f, calLog, calG);
        end
end
r.cfg = c; r.seed = seed; r.time = toc(t0);
r.score = cpr_score(X, xh, ns);
r.scoreVec = vec_score(X, xh, c.nd, c.npol);
info.snrCombinedDb = 10*log10(mean(abs(z).^2) / wvar);   % realized SNR of the combined signals (round 7)
r.info = info; r.R = S.R; r.m = S.m; r.n = S.n;
if isempty(xh), r.Xh = []; else, r.Xh = reshape(xh, size(X)); end

    function calReset()
        calD = zeros(S.ent.nPorts, 1); calG = calD; calCalls = 0; calLog = [];
        A.setPhi(PhiSolver);
        reset0();
    end

    function est = calHook(est)
        calCalls = calCalls + 1;
        if calCalls >= c.selfCalStart
            sigL = c.selfCalPriorDeg / 360 * (299792458 / mean(f) * 100);   % cm at the band centre
            calOpts.g = calG;
            [Phi1, calD, ci, calG] = cpr_selfcal(S.ent, A.syn(est.xhat), y, S.lam, calD, c.selfCalGN, sigL, calOpts);
            A.setPhi(Phi1);
            calLog = [calLog; ci.res0_dB, ci.res_dB]; %#ok<AGROW>
        end
        if ~isempty(gauge0), est = gauge0(est); end
    end
end

function s = cal_score(d, S, f, calLog, g)
% error of the estimated length corrections (cm, and rms phase at the band centre in degrees), with the
% common length removed from both; gain: rms error of the log-amplitudes in dB, common gain removed
s.d = d; s.log = calLog;
s.g = g; s.rmsGainTrueDb = NaN; s.rmsGainErrDb = NaN;
if isfield(S, 'gTrue') && ~isempty(S.gTrue) && any(S.gTrue)
    gt = S.gTrue - mean(S.gTrue); s.rmsGainTrueDb = 20/log(10) * rms(gt);
    s.rmsGainErrDb = 20/log(10) * rms(gt - (g - mean(g)));
end
if isempty(S.dLtrue), s.rmsErrDeg = NaN; s.rmsTrueDeg = NaN; return, end
lamC = 299792458 / mean(f) * 100;
t = S.dLtrue - mean(S.dLtrue);
s.rmsTrueDeg = 360 * rms(t) / lamC;
s.rmsErrDeg = 360 * rms(t - (d - mean(d))) / lamC;
end

function s = vec_score(X, xh, nd, npol)
% NMSE of the full (all-polarization) pattern: one phase per wavelength shared by the
% polarizations, and one global phase.
if isempty(xh), s = struct('raw', NaN, 'perlambda', NaN, 'global', NaN, 'each', NaN(1, nd)); return, end
X = reshape(X, [], nd, npol); Xh = reshape(xh, [], nd, npol);
nr = @(E) 20*log10(norm(E(:))/norm(X(:)));
s.raw = nr(X - Xh);
s.global = nr(X - Xh*exp(1i*angle(Xh(:)'*X(:))));
s.each = zeros(1, nd);
for k = 1:nd
    a = Xh(:,k,:); b = X(:,k,:);
    Xh(:,k,:) = a * exp(1i*angle(a(:)'*b(:)));
    d = Xh(:,k,:) - b; s.each(k) = 20*log10(norm(d(:))/norm(b(:)));
end
s.perlambda = nr(X - Xh);
end
