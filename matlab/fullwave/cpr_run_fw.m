function r = cpr_run_fw(c, seed, D, legMode)
%CPR_RUN_FW  One trial on the Section VI full-wave IC geometry.
%   r = cpr_run_fw(cpr_cfg('meas','complex','P',4,'K',4,'ratio',0.15,'nlegs',3,'snrdB',50), 1, D, 'real')
% D = fw_data(...) (load once, reuse). Same solvers, basis, gauge and scoring as
% cpr_run.m; the pattern grid is 16 theta x 32 rotations x nf (D sizes, not cfg nw/nd).
cpr_setup();                                       % puts GAMPmatlab (main, phase) on the path
if nargin < 4, legMode = 'real'; end
[nTh, nPh, nf] = size(D.Xring);
S = make_readings_fw(c, D, seed, legMode);          % calls rng(seed)
z = S.z;
wvar = cpr_noisevar(c, z, D.Xring);
w = sqrt(wvar/2) * (randn(S.m, 2) * [1; 1i]);
A = SparseDftLinTrans(S.Phi, nTh, nPh, nf, cpr_basis(c.basis, D.f, 0), c.sqMode);
t0 = tic;
switch c.meas
    case 'complex'
        o = struct();
        if c.gauge && nf > 1, o.gauge = @(e) cpr_gauge(e, A, 'delay', D.f, c.basis, S.rowTone); end
        [xh, info] = cpr_gamp_cs(z + w, A, o);
    case 'phaseless'
        o = struct('minTry', c.minTry, 'maxTry', c.maxTry, 'nitEM', c.nitEM, 'SNRdB', min(c.solverSNRdB, c.snrdB));
        o.nitGAMP = cpr_settings().fullwave.nitGAMP;           % 200 in the paper's Section VI runs
        o.reset = @() A.setM(cpr_basis(c.basis, D.f, 0));
        if c.gauge && nf > 1
            if c.K == 1, mode = 'perlambda'; else, mode = 'delay'; end
            o.gauge = @(e) cpr_gauge(e, A, mode, D.f, c.basis, S.rowTone);
        end
        [xh, info] = feval(c.solver, abs(z + w), A, o);
    otherwise
        error('cpr_run_fw: unknown meas %s', c.meas);
end
r.cfg = c; r.seed = seed; r.legMode = legMode; r.time = toc(t0);
r.score = cpr_score(D.Xring, xh, nf);
r.info = info; r.R = S.R; r.m = S.m; r.n = S.n; r.legMismatch_dB = S.legMismatch_dB;
if isempty(xh), r.Xh = single([]); else, r.Xh = single(reshape(xh, nTh, nPh, nf)); end
end
