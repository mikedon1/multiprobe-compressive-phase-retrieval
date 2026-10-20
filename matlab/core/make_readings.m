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
%        rotations; nlegs leg probes with their own cables. c.legDir = 'phi' (default):
%        leg q sits at the anchor row theta = 90 deg (row nw/4+1) of the rotation q steps
%        ahead, wrapping (as in the full-wave model); 'theta': legs at theta rows
%        nw-nlegs+1..nw of the NEXT rotation (earlier versions).
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
        legDir = 'theta'; if isfield(c, 'legDir'), legDir = c.legDir; end
        if strcmp(legDir, 'phi')
            if nl >= nh, error('make_readings: %d phi legs need more than %d rotations', nl, nh); end
            legRow = repmat(nw/4 + 1, nl, 1); legOff = (1:nl).';   % anchor theta = 90 deg, q rotations ahead
        else
            legRow = (nw - nl + 1 : nw).'; legOff = ones(nl, 1);   % theta column of the next rotation
        end
        % IC with a theta-stepped ring for calibration (2026-10-19): the last fraction icShiftFrac of the readings of
        % every rotation are taken with the ring turned by icShift theta samples (an arch offset, as in StarLab), so
        % each port sees two directions; default 0 = plain IC
        % c.icShiftMode 'perRot' (default): the last fraction of EVERY rotation's readings (ring out and back at each
        % phi); 'end': the last fraction of ALL readings (scan normally, turn the ring once, finish the scan)
        sh = [0 0 0];
        if isfield(c, 'icShiftFrac') && c.icShiftFrac > 0
            if isfield(c, 'icShiftMode') && strcmp(c.icShiftMode, 'end'), sh = [c.icShift, round(c.icShiftFrac*R), 1];
            else, sh = [c.icShift, round(c.icShiftFrac*perRot), 0]; end
        end
        angOf = @(r, idx) ring_leg_ang(r, idx, c.arch, pos0, ntx, legRow, legOff, perRot, nw, nh, sh);
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
EP = I; ET = I;                                           % port and tone index of every entry (self-calibration)
rowTone = zeros(mTot, 1);
ptr = 0; row = 0;
gauss = isfield(c, 'weights') && strcmp(c.weights, 'gauss');
% combiner normalization (round 7, 2026-10-20): 'unit' (default) = unit branch transfers, an idealized buffered/lossless
% combiner; 'passive' = ideal passive matched P-way combiner of the selected ports, branch amplitude 1/sqrt(P)
passiveComb = isfield(c, 'combNorm') && strcmp(c.combNorm, 'passive');
if gauss, gs = RandStream('mt19937ar', 'Seed', seed + 7919); end   % separate stream: other draws unchanged
% calibration error: the solver assumes the calibrated cable phases, the data use the true ones.
% c.calErrModel = 'length' (default, physical): each port's cable length is off by a fixed dL_p, so its
% phase error 2*pi*dL_p/lambda_k is fixed per port and grows with frequency; dL_p ~ N(0, sigma_L^2) with
% sigma_L set so that the rms phase error at the band centre is calErrDeg. 'tone' (earlier versions):
% independent N(0, calErrDeg^2) per port and tone. Separate stream, so all other draws are unchanged.
% 'drift' (2026-10-19): a thermal drift, every cable longer by the same FRACTION, dL_p = a*L_p, with a scaled so
% that the rms phase error over the ports at the band centre is calErrDeg (random sign).
% c.calGainDb > 0 (2026-10-19): also a fixed gain error per port, N(0, calGainDb^2) in dB, the same at every tone
% (cable/switch loss), from its own stream so the length draws are unchanged.
phDeg = 0; if isfield(c, 'calErrDeg'), phDeg = c.calErrDeg; end
gDb = 0; if isfield(c, 'calGainDb'), gDb = c.calGainDb; end
calErr = phDeg > 0 || gDb > 0;
if calErr
    es = RandStream('mt19937ar', 'Seed', seed + 104729);
    mdl = 'length'; if isfield(c, 'calErrModel'), mdl = c.calErrModel; end
    lamC = 299792458 / mean(f) * 100;                                         % band-centre wavelength (cm)
    switch mdl
        case 'length'
            dL = (phDeg/360) * lamC * randn(es, nPool*npol, 1);               % cm, one per port
            E = exp(1i*2*pi*dL ./ lam);                                       % true / assumed, port x tone
        case 'drift'
            Lp = L(:); a = sign(randn(es)) * (phDeg/360) * lamC / rms(Lp);    % fractional length change
            dL = a * Lp;
            E = exp(1i*2*pi*dL ./ lam);
        case 'tone'
            E = exp(1i*(phDeg*pi/180)*randn(es, nPool*npol, nd));             % true / assumed, per port and tone
        otherwise
            error('make_readings: unknown calErrModel %s', mdl);
    end
    gTrue = zeros(nPool*npol, 1);
    if gDb > 0
        gsr = RandStream('mt19937ar', 'Seed', seed + 104743);
        gTrue = (gDb/20) * log(10) * randn(gsr, nPool*npol, 1);              % natural-log amplitude per port
        E = E .* exp(gTrue);
    end
end
% timing jitter: the detector samples at t + dt (dt ~ N(0, jitterPs^2), independent per sample) while the
% solver assumes t. Separate stream, so all other draws are unchanged.
jit = phaseless && isfield(c, 'jitterPs') && c.jitterPs > 0;
if jit, js = RandStream('mt19937ar', 'Seed', seed + 15485863); end
mism = calErr || jit;                                     % data and solver matrices differ
% tone window (2026-10-19): each reading draws its K tones from a window of span adjacent frequency bins at a
% random position, so the detector only sees beats up to (span - 1)*df. span >= nd (default) = the whole band.
span = nd; if isfield(c, 'toneSpan'), span = min(c.toneSpan, nd); end
if span < nd && ~phaseless, error('make_readings: toneSpan is implemented for phaseless readings only'); end
if mism, VA = complex(V); end
for r = 1:R
    idx = randperm(nPool*npol, P).';                      % active ports
    pr = mod(idx - 1, nPool) + 1; q = floor((idx - 1)/nPool) + 1;   % probe, polarization
    a = angOf(r, pr);                                     % P x 1 grid points
    if span < nd                                          % tone window (span adjacent bins, may overhang the
        s0 = randi(span + nd - 1) - span + 1;             %   band edges so every frequency is covered equally;
        win = max(1, s0):min(nd, s0 + span - 1);          %   out-of-band tones are dropped)
        tones = sort(win(randperm(numel(win), min(K, numel(win)))));
    else
        tones = sort(randperm(nd, K));
    end
    nt = numel(tones);                                    % tones in this reading (K except near the edges)
    W = exp(1i*2*pi*L(sub2ind(size(L), pr, q)) ./ lam(tones));   % P x nt cable phases
    if gauss, W = W .* ((randn(gs, P, 1) + 1i*randn(gs, P, 1)) / sqrt(2)); end   % CN(0,1), same for all tones
    if passiveComb, W = W / sqrt(P); end                  % ideal passive P-way combiner (known; data and model)
    cols = a + (tones - 1)*n2 + (q - 1)*nd*n2;            % P x K pattern columns
    if calErr, Et = E(idx, tones); else, Et = ones(P, nt); end   % true/assumed phase, P x nt
    if ~phaseless
        for q = 1:K
            row = row + 1;
            I(ptr+(1:P)) = row; J(ptr+(1:P)) = cols(:,q); V(ptr+(1:P)) = W(:,q) .* Et(:,q);
            EP(ptr+(1:P)) = idx; ET(ptr+(1:P)) = tones(q);
            if mism, VA(ptr+(1:P)) = W(:,q); end
            ptr = ptr + P; rowTone(row) = tones(q);
        end
    else
        psi = 2*pi*rand(1, nt);
        t = rand(K, 1) / df;                              % K time samples in one envelope period
        tT = t; if jit, tT = t + c.jitterPs*1e-12*randn(js, K, 1); end   % true sample times
        for q = 1:K
            row = row + 1;
            ph = exp(1i*(2*pi*(f(tones) - f(1))*t(q) + psi));   % 1 x nt (assumed)
            Wq = W .* ph;
            WqT = Wq .* Et;
            if jit, WqT = W .* exp(1i*(2*pi*(f(tones) - f(1))*tT(q) + psi)) .* Et; end
            I(ptr+(1:P*nt)) = row; J(ptr+(1:P*nt)) = cols(:); V(ptr+(1:P*nt)) = WqT(:);
            EP(ptr+(1:P*nt)) = repmat(idx, nt, 1); ET(ptr+(1:P*nt)) = reshape(repmat(tones, P, 1), [], 1);
            if mism, VA(ptr+(1:P*nt)) = Wq(:); end
            ptr = ptr + P*nt;
            if K == 1, rowTone(row) = tones; end
        end
    end
end
S.Phi = sparse(I(1:ptr), J(1:ptr), V(1:ptr), mTot, n);                    % data (true phases)
if mism, S.PhiA = sparse(I(1:ptr), J(1:ptr), VA(1:ptr), mTot, n); end   % solver (calibrated phases, nominal times)
S.rowTone = rowTone; S.R = R; S.m = mTot; S.n = n;
S.f = f; S.lam = lam; S.df = df; S.L = L;
% entries of the solver's matrix with their port and tone (cpr_selfcal re-phases them per port)
vA = V(1:ptr); if mism, vA = VA(1:ptr); end
S.ent = struct('i', I(1:ptr), 'j', J(1:ptr), 'v', vA, 'port', EP(1:ptr), 'tone', ET(1:ptr), ...
    'nPorts', nPool*npol, 'm', mTot, 'n', n);
S.dLtrue = []; if calErr && any(strcmp(mdl, {'length', 'drift'})), S.dLtrue = dL; end   % true length errors (cm), for scoring
S.gTrue = []; if calErr, S.gTrue = gTrue; end            % true log-amplitude errors per port (0 if no gain error)
end

function a = ring_leg_ang(r, idx, arch, pos0, ntx, legRow, legOff, perRot, nw, nh, sh)
% grid points seen by ring/leg probes idx at reading r (1-based); sh = [theta step, readings per rotation with
% the stepped ring] (IC only)
q = floor((r-1)/perRot);                                  % phi rotation 0..nh-1
a = zeros(numel(idx), 1);
isRing = idx <= ntx;
ir = idx(isRing);
if strcmp(arch, 'IC')
    th = pos0(ir);
    if sh(2) > 0
        if sh(3) == 1, turned = r > perRot*nh - sh(2); else, turned = mod(r-1, perRot) >= perRot - sh(2); end
        if turned, th = mod(th - 1 + sh(1), nw) + 1; end
    end
else
    th = mod(pos0(ir) - 1 - (r-1), nw) + 1;               % ring steps one theta sample per reading
end
a(isRing) = th + q*nw;
il = idx(~isRing) - ntx;
a(~isRing) = legRow(il) + mod(q + legOff(il), nh)*nw;    % legOff rotations ahead, wrapping
end
