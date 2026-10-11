function r = cpr_run_pattern(c, seed, X, f)
%CPR_RUN_PATTERN  cpr_run.m for an arbitrary pattern X (nw x nh x nd*npol) and frequencies f.
%Used by the IFT comparison (near-field data). c.nw and c.nh are set from X, so the grid need
%not be nh = nw/2. With c.npol = 2 (dual polarization, Section V), X stacks the polarizations
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
wvar = 10^(-c.snrdB/10) * mean(abs(z).^2);
w = sqrt(wvar/2) * (randn(S.m, 2) * [1; 1i]);
A = SparseDftLinTrans(S.Phi, c.nw, c.nh, ns, cpr_basis(c.basis, f, 0, c.npol), c.sqMode);
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
        o.reset = @() A.setM(cpr_basis(c.basis, f, 0, c.npol));
        if c.gauge && c.nd > 1
            if c.K == 1, mode = 'perlambda'; else, mode = 'delay'; end
            o.gauge = @(e) cpr_gauge(e, A, mode, f, c.basis, S.rowTone);
        end
        [xh, info] = feval(c.solver, abs(z + w), A, o);
end
r.cfg = c; r.seed = seed; r.time = toc(t0);
r.score = cpr_score(X, xh, ns);
r.scoreVec = vec_score(X, xh, c.nd, c.npol);
r.info = info; r.R = S.R; r.m = S.m; r.n = S.n;
if isempty(xh), r.Xh = []; else, r.Xh = reshape(xh, size(X)); end
end

function s = vec_score(X, xh, nd, npol)
% NMSE of the full (all-polarization) pattern: one phase per wavelength shared by the
% polarizations, and one global phase.
if isempty(xh), s = struct('perlambda', NaN, 'global', NaN, 'each', NaN(1, nd)); return, end
X = reshape(X, [], nd, npol); Xh = reshape(xh, [], nd, npol);
nr = @(E) 20*log10(norm(E(:))/norm(X(:)));
s.global = nr(X - Xh*exp(1i*angle(Xh(:)'*X(:))));
s.each = zeros(1, nd);
for k = 1:nd
    a = Xh(:,k,:); b = X(:,k,:);
    Xh(:,k,:) = a * exp(1i*angle(a(:)'*b(:)));
    d = Xh(:,k,:) - b; s.each(k) = 20*log10(norm(d(:))/norm(b(:)));
end
s.perlambda = nr(X - Xh);
end
