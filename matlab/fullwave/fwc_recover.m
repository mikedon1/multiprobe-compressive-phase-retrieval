function r = fwc_recover(c, seed, D)
%FWC_RECOVER  One Section VI trial with embedded-probe calibration: recover the AUT's FREE-SPACE mode
%coefficients from phaseless readings of the dual-polarized ring (no legs).
%   r = fwc_recover(cpr_cfg('arch', 'IC', 'ratio', 0.6, 'snrdB', 50), seed, D)
% D fields:
%   Xring  np x nrot x nf  measured AUT->probe transmissions (HFSS full model), ports in map order
%   R      np x nmode x nf  calibrated embedded responses of the same ports (fwc_probe_resp)
%   modes  nmode x 3 [l m type];  psi (rad, 1 x nrot);  f (Hz, 1 x nf)
%   aTrue  nmode x nf  free-space mode coefficients of the AUT (Solve A, fwc_aut_modes .a)
% Readings: make_readings_fw (on/off switches, P ports, K tones, U(0, cableMaxCm) cables, no legs), data =
% Phi * Xring(:) + noise (the real full-wave responses). Model: CalVshLinTrans (Phi, G_f = R_f exp(-i m psi),
% DCT + delay along frequency), PR-GAMP with the default structured prior and delay re-expression.
% Scores (dB, full sphere: the VSH are orthonormal over the sphere, so the coefficient error is the field
% error): r.nmse.perfreqScale removes one complex factor per frequency (relative pattern at each frequency,
% the main score); r.nmse.perfreqPhase removes one phase per frequency; r.nmse.globalPhase / globalScale one phase /
% one complex factor for all frequencies (tests the phase between frequencies); r.perfreqFactor the fitted factors; r.nmse.model is the data-domain
% model error of the truth, || G aTrue - Xring || / || Xring || with one complex factor per frequency.
cpr_setup();
[np, nrot, nf] = size(D.Xring); nm = size(D.modes, 1);
G = cell(1, nf);
for f = 1:nf
    E = exp(-1i * D.modes(:, 2).' .* D.psi(:));                     % nrot x nmode
    Gf = zeros(np*nrot, nm);
    for l = 1:nrot
        Gf((l-1)*np + (1:np), :) = D.R(:, :, f) .* E(l, :);
    end
    G{f} = Gf;
end
Dm = struct('Xring', D.Xring, 'f', D.f(:).', 'nlegsAvail', 0, 'Xleg', [], 'anchor', 1, 'legShift', 1);
c.nlegs = 0;
S = make_readings_fw(c, Dm, seed, 'ideal');
z = S.Phi * D.Xring(:);
wvar = cpr_noisevar(c, z, D.Xring);
w = sqrt(wvar/2) * (randn(S.m, 2) * [1; 1i]);
A = CalVshLinTrans(S.Phi, G, cpr_basis(c.basis, D.f, 0));
o = struct('minTry', c.minTry, 'maxTry', c.maxTry, 'nitEM', c.nitEM, 'SNRdB', min(c.solverSNRdB, c.snrdB), ...
    'nitGAMP', c.nitGAMP, 'groupPrior', c.groupPrior);
if isfield(c, 'verbose'), o.verbose = c.verbose; end
o.reset = @() A.setM(cpr_basis(c.basis, D.f, 0));
if c.gauge && nf > 1, o.gauge = @(e) cpr_gauge(e, A, 'delay', D.f, c.basis, S.rowTone); end
t0 = tic;
[ah, info] = feval(c.solver, abs(z + w), A, o);
r.cfg = c; r.seed = seed; r.time = toc(t0); r.info = info; r.m = S.m; r.n = numel(D.Xring);
ah = reshape(ah, nm, nf); r.aHat = ah;
r.nmse.perfreqScale = score(D.aTrue, ah, 'scale');
r.nmse.perfreqPhase = score(D.aTrue, ah, 'phase');
r.nmse.globalPhase = score(D.aTrue(:), ah(:), 'phase');   % one phase for all frequencies (2026-10-19)
r.nmse.globalScale = score(D.aTrue(:), ah(:), 'scale');   % one complex factor for all frequencies
[~, r.perfreqFactor] = score(D.aTrue, ah, 'scale');      % fitted factor per frequency (smooth in f? Sec. VII)
Xm = reshape(A.probes(D.aTrue(:)), np*nrot, nf); Xd = reshape(D.Xring, np*nrot, nf);
r.nmse.model = score(Xd, Xm, 'scale');
end

function [e, kf] = score(T, Y, mode)
% NMSE (dB) of Y against T (columns = frequencies) after one complex factor ('scale') or one phase
% ('phase') per column
err = 0; kf = zeros(1, size(T, 2));
for f = 1:size(T, 2)
    t = T(:, f); y = Y(:, f); k = (y' * t) / max(y' * y, eps);
    if strcmp(mode, 'phase'), k = exp(1i * angle(k)); end
    kf(f) = k;
    err = err + norm(k * y - t)^2;
end
e = 10*log10(err / norm(T, 'fro')^2);
end
