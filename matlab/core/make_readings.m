function S = make_readings(c, f, seed)
%MAKE_READINGS  Pattern-domain sensing matrix for one trial.
%
%   S = make_readings(c, f, seed)   c = cpr_cfg(...), f = tone frequencies (Hz, 1 x nd)
%
% Physical model: on/off switches. One READING = one switch setting (P
% active probes) while the AUT transmits K tones drawn at random from the nd
% frequencies. Probe p at tone k contributes exp(1i*2*pi*L_p/lambda_k)
% (cable length L_p ~ U(0, cableMaxCm) cm, fixed to the physical probe and
% carried with it through rotations).
%   complex  : the coherent receiver FFTs the record -> K rows per reading,
%              one per tone, all with the same probe set.
%   phaseless: a fast square-law detector is sampled at T = K random times t
%              in one envelope period 1/df -> K rows per reading, each
%              sum_k exp(1i*(2*pi*(f_k - f_1)*t + psi_k)) * (tone-k weights),
%              with known random tone phases psi_k. y = |row*x + noise|.
%              K = 1 is the paper's original single-frequency model.
% Architectures (pattern grid nw x nh, theta 0-360 x phi 0-180):
%   IS : a probe at every grid point; no rotation.
%   RS : round(rnRS*nh) probes per theta row at random phi positions; the
%        whole set steps one phi sample per reading (RS convention).
%   IC : a full theta ring (nw probes); readings grouped into nh phi
%        rotations; nlegs leg probes at theta rows nw-nlegs+1..nw sample the
%        NEXT phi rotation (wrapping to the first), with their own cables.
%   RC : round(rnRC*nw) ring probes at random theta positions; the ring steps
%        one theta sample per reading and readings are grouped into nh phi
%        rotations (RC convention); legs as IC.
% Active probes are drawn from the whole installed pool (ring + legs).
% Dual polarization (c.npol = 2): every installed probe has two ports
% (theta and phi), each with its own calibrated cable, both on the switch network. A reading
% closes P PORTS drawn from all ports, so one reading can combine the two polarizations; the
% cross terms then measure the phase between them (as multi-tone readings do between
% wavelengths). The pattern stacks the polarizations: slices [pol 1: k = 1..nd, pol 2: ...].
% npol = 1 reproduces the single-polarization generator exactly (same random draws).
%
% Columns are pattern samples a + (k-1)*n2 + (q-1)*nd*n2, a = itheta + (iphi-1)*nw, q = pol.
% Output S: Phi (sparse m x n), rowTone (tone index per row for complex or
% K=1, else 0), R (readings), m, n, f, lam (cm), df.
rng(seed);
npol = 1; if isfield(c, 'npol'), npol = c.npol; end
nw = c.nw; nh = c.nh; nd = c.nd; n2 = nw*nh; n = n2*nd*npol;
K = min(c.K, nd); P = c.P;
f = f(:).'; lam = 299792458 ./ f * 100;
if nd > 1, df = f(2) - f(1); else, df = 50e6; end
m0 = round(c.ratio * n);
R = max(1, round(m0 / K));
if any(strcmp(c.arch, {'IC', 'RC'})), R = max(nh, round(R/nh)*nh); end

% ---------- installed probes: cable length and sampled grid point per reading
switch c.arch
    case 'IS'
        L = rand(n2, 1) * c.cableMaxCm;
        nPool = n2;
        angOf = @(r, idx) idx(:);                         % probe index = grid point
    case 'RS'
        npr = round(c.rnRS * nh);
        if npr < 1, error('make_readings: rnRS*nh < 1'); end
        th = zeros(nw*npr, 1); ph0 = th;
        for i = 1:nw
            j = randperm(nh, npr);
            th((i-1)*npr + (1:npr)) = i; ph0((i-1)*npr + (1:npr)) = j;
        end
        L = rand(nw*npr, 1) * c.cableMaxCm;
        nPool = nw*npr;
        angOf = @(r, idx) th(idx) + (mod(ph0(idx) - 1 + (r-1), nh)) * nw;
    case {'IC', 'RC'}
        nl = c.nlegs;
        if strcmp(c.arch, 'IC')
            ntx = nw; pos0 = (1:nw).';
        else
            ntx = round(c.rnRC * nw); pos0 = randperm(nw, ntx).';
        end
        L = rand(ntx + nl, 1) * c.cableMaxCm;            % ring then legs
        nPool = ntx + nl;
        perRot = R / nh;
        legRow = (nw - nl + 1 : nw).';
        angOf = @(r, idx) ring_leg_ang(r, idx, c.arch, pos0, ntx, legRow, perRot, nw, nh);
    otherwise
        error('make_readings: unknown arch %s', c.arch);
end
if npol > 1, L = [L, rand(nPool, npol - 1) * c.cableMaxCm]; end   % one cable per port
if P > nPool*npol, error('make_readings: P=%d exceeds the %d installed ports', P, nPool*npol); end

% ---------- readings -> sparse rows
phaseless = strcmp(c.meas, 'phaseless');
nnzRow = P*K*phaseless + P*(~phaseless);
mTot = R*K;
I = zeros(mTot*nnzRow, 1); J = I; V = complex(I);
rowTone = zeros(mTot, 1);
ptr = 0; row = 0;
gauss = isfield(c, 'weights') && strcmp(c.weights, 'gauss');
if gauss, gs = RandStream('mt19937ar', 'Seed', seed + 7919); end   % separate stream: other draws unchanged
% calibration error: the solver assumes the calibrated cable phases, the data use phases
% that differ by N(0, calErrDeg^2) per port and tone. Separate stream, so all other draws are unchanged.
calErr = isfield(c, 'calErrDeg') && c.calErrDeg > 0;
if calErr
    es = RandStream('mt19937ar', 'Seed', seed + 104729);
    E = exp(1i*(c.calErrDeg*pi/180)*randn(es, nPool*npol, nd));   % true / assumed, per port and tone
    VA = complex(V);
end
for r = 1:R
    idx = randperm(nPool*npol, P).';                      % active ports
    pr = mod(idx - 1, nPool) + 1; q = floor((idx - 1)/nPool) + 1;   % probe, polarization
    a = angOf(r, pr);                                     % P x 1 grid points
    tones = sort(randperm(nd, K));
    W = exp(1i*2*pi*L(sub2ind(size(L), pr, q)) ./ lam(tones));   % P x K cable phases
    if gauss, W = W .* ((randn(gs, P, 1) + 1i*randn(gs, P, 1)) / sqrt(2)); end   % CN(0,1), same for all tones
    cols = a + (tones - 1)*n2 + (q - 1)*nd*n2;            % P x K pattern columns
    if calErr, Et = E(idx, tones); else, Et = ones(P, K); end   % true/assumed phase, P x K
    if ~phaseless
        for q = 1:K
            row = row + 1;
            I(ptr+(1:P)) = row; J(ptr+(1:P)) = cols(:,q); V(ptr+(1:P)) = W(:,q) .* Et(:,q);
            if calErr, VA(ptr+(1:P)) = W(:,q); end
            ptr = ptr + P; rowTone(row) = tones(q);
        end
    else
        psi = 2*pi*rand(1, K);
        t = rand(K, 1) / df;                              % K time samples in one envelope period
        for q = 1:K
            row = row + 1;
            ph = exp(1i*(2*pi*(f(tones) - f(1))*t(q) + psi));   % 1 x K
            Wq = W .* ph;
            WqT = Wq .* Et;
            I(ptr+(1:P*K)) = row; J(ptr+(1:P*K)) = cols(:); V(ptr+(1:P*K)) = WqT(:);
            if calErr, VA(ptr+(1:P*K)) = Wq(:); end
            ptr = ptr + P*K;
            if K == 1, rowTone(row) = tones; end
        end
    end
end
S.Phi = sparse(I(1:ptr), J(1:ptr), V(1:ptr), mTot, n);                    % data (true phases)
if calErr, S.PhiA = sparse(I(1:ptr), J(1:ptr), VA(1:ptr), mTot, n); end   % solver (calibrated phases)
S.rowTone = rowTone; S.R = R; S.m = mTot; S.n = n;
S.f = f; S.lam = lam; S.df = df; S.L = L;
end

function a = ring_leg_ang(r, idx, arch, pos0, ntx, legRow, perRot, nw, nh)
% grid points seen by ring/leg probes idx at reading r (1-based)
q = floor((r-1)/perRot);                                  % phi rotation 0..nh-1
a = zeros(numel(idx), 1);
isRing = idx <= ntx;
ir = idx(isRing);
if strcmp(arch, 'IC')
    th = pos0(ir);
else
    th = mod(pos0(ir) - 1 - (r-1), nw) + 1;               % ring steps one theta sample per reading
end
a(isRing) = th + q*nw;
il = idx(~isRing) - ntx;
a(~isRing) = legRow(il) + mod(q + 1, nh)*nw;              % next phi rotation, wrapping
end
