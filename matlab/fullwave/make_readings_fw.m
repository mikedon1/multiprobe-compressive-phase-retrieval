function S = make_readings_fw(c, D, seed, legMode)
%MAKE_READINGS_FW  Readings for the Section VI full-wave IC geometry.
%
%   S = make_readings_fw(c, D, seed, legMode)
%     c       cpr_cfg(...) -- uses meas, P, K, ratio, nlegs (<= D.nlegsAvail),
%             cableMaxCm
%     D       struct with the fields listed below (fwc_run_trials)
%     legMode 'real'  : data use each leg's OWN response (real ports / exact
%                       position); the recovery model maps leg k at rotation L
%                       to the anchor at rotation L + D.legShift*k. Any leg
%                       coupling is therefore model mismatch, as in a real
%                       system that never knows the leg's own "true" values.
%             'ideal' : data built from that same mapping (no mismatch).
%
% Same physical model as make_readings.m: on/off switches, P active probes from
% the pool (16 ring + nlegs legs), K random tones per reading, calibrated cable
% lengths U(0, cableMaxCm) carried by each physical probe, K per-tone values
% (complex) or K time samples of the square-law detector (phaseless). Readings
% are split evenly over the 32 rotations.
%
% Output: S.Phi   model matrix (m x n, n = 16*32*nf) acting on Xring(:)
%         S.z     noiseless data (m x 1) -- includes the real-leg mismatch
%         S.rowTone, S.R, S.m, S.n, S.f, S.lam, S.df
rng(seed);
[nTh, nPh, nf] = size(D.Xring); n2 = nTh*nPh; n = n2*nf;
nl = min(c.nlegs, D.nlegsAvail);
if nl == 0 && c.nlegs > 0 && isempty(D.Xleg)
    nl = min(c.nlegs, 3);                       % 17-port data: legs only from the mapping
    if strcmp(legMode, 'real'), error('make_readings_fw:noLegs', 'D has no leg ports; use legMode=''ideal'''); end
end
K = min(c.K, nf); P = c.P;
f = D.f; lam = 29979245800 ./ f;
if nf > 1, df = f(2) - f(1); else, df = 50e6; end
R = max(nPh, round(round(c.ratio*n/K)/nPh)*nPh); perRot = R/nPh;
L = rand(nTh + nl, 1) * c.cableMaxCm;            % ring cables then leg cables
nPool = nTh + nl;
if P > nPool, error('make_readings_fw:P', 'P=%d exceeds %d probes', P, nPool); end
phaseless = strcmp(c.meas, 'phaseless');
nnzRow = P*(phaseless*K + ~phaseless);
mTot = R*K;
% Leg calibration in the OPERATOR: if D.legCal (nl x nf) is given, the
% recovery model's leg entries carry the calibration factor, w * kappa_q(f) * x_anchor, while the
% simulated (physical, switch-only) data keep each leg's actual response w * l_q. Nothing is
% divided before combination.
hasCal = isfield(D, 'legCal') && ~isempty(D.legCal);
I = zeros(mTot*nnzRow, 1); J = I; V = complex(I); Vd = V; isLeg = false(size(I)); legReal = complex(I);
rowTone = zeros(mTot, 1); ptr = 0; row = 0;
for r = 1:R
    q = floor((r-1)/perRot);                          % rotation 0..nPh-1
    idx = randperm(nPool, P).';
    tones = sort(randperm(nf, K));
    a = zeros(P, 1); legk = zeros(P, 1);
    ring = idx <= nTh;
    a(ring) = idx(ring) + q*nTh;
    legk(~ring) = idx(~ring) - nTh;
    a(~ring) = D.anchor + mod(q + D.legShift*legk(~ring), nPh)*nTh;   % model column for each leg
    W = exp(1i*2*pi*L(idx) ./ lam(tones));            % P x K, physical (data) weights
    Wm = W;                                           % model weights
    if hasCal && any(~ring)
        Wm(~ring, :) = W(~ring, :) .* D.legCal(legk(~ring), tones);
    end
    cols = a + (tones - 1)*n2;
    if strcmp(legMode, 'real') && any(~ring)          % each leg's own response, per tone
        lr = zeros(P, K);
        for kk = 1:K
            lr(~ring, kk) = D.Xleg(legk(~ring), q+1, tones(kk));
        end
    else
        lr = zeros(P, K);
    end
    if ~phaseless
        for kk = 1:K
            row = row + 1; s = ptr + (1:P);
            I(s) = row; J(s) = cols(:,kk); V(s) = Wm(:,kk); Vd(s) = W(:,kk);
            isLeg(s) = ~ring; legReal(s) = lr(:,kk);
            ptr = ptr + P; rowTone(row) = tones(kk);
        end
    else
        psi = 2*pi*rand(1, K); t = rand(K, 1) / df;
        for qq = 1:K
            row = row + 1; s = ptr + (1:P*K);
            ph = exp(1i*(2*pi*(f(tones) - f(1))*t(qq) + psi));
            Wq = Wm .* ph; Wqd = W .* ph;
            I(s) = row; J(s) = cols(:); V(s) = Wq(:); Vd(s) = Wqd(:);
            isLeg(s) = repmat(~ring, K, 1); legReal(s) = lr(:);
            ptr = ptr + P*K;
            if K == 1, rowTone(row) = tones; end
        end
    end
end
I = I(1:ptr); J = J(1:ptr); V = V(1:ptr); Vd = Vd(1:ptr); isLeg = isLeg(1:ptr); legReal = legReal(1:ptr);
S.Phi = sparse(I, J, V, mTot, n);
x = D.Xring(:);
z = S.Phi * x;
if strcmp(legMode, 'real')
    % physical data: each leg entry contributes w * (its own response); the model has (calibrated) w * x_anchor
    dz = accumarray(I(isLeg), Vd(isLeg) .* legReal(isLeg) - V(isLeg) .* x(J(isLeg)), [mTot 1]);
    z = z + dz;
    S.legMismatch_dB = 20*log10(norm(dz)/norm(z));
else
    S.legMismatch_dB = -Inf;
end
S.z = z; S.rowTone = rowTone; S.R = R; S.m = mTot; S.n = n;
S.f = f; S.lam = lam; S.df = df; S.L = L; S.nlegs = nl;
end
