function T = jac_rot_geometry(seed, tols)
%JAC_ROT_GEOMETRY  (2026-10-25, study 3 of paper/round10/MORE_SIMS_PLAN.md) Local amplitude-Jacobian sensitivity to rotation phases,
%compared across sampling geometries of the Section VII pipeline (ideal far-field dual-polarized point probes, dipole AUT coefficients,
%base readings: IC, 16 frequencies, K = 4, m/n = 0.3, L = 12, one seed). For every geometry:
%  Js   real Jacobian of the amplitudes |z_i| with respect to the real and imaginary parts of the free-space mode coefficients on the
%       support that holds a fraction tol of the energy (the true support is used, so this shows what sparsity can contribute, not what
%       a solver finds)
%  JG   Jacobian with respect to rotation phases g_l (the AUT data of rotation l multiplied by exp(i g_l)): (a) one phase per rotation,
%       (b) one phase per rotation and frequency. The all-ones direction (the global phase) is removed.
%  margin = singular values of (I - P_Js) JG Q normalized by JG Q (as round7_checks part 4c for cable lengths): the smallest is the
%       rotation-phase direction that the amplitudes can least distinguish from a change of the coefficients, in dB (0 dB = fully
%       visible, -inf = invisible). count = directions below -30 dB.
% LOCAL, support-conditioned evidence only: no uniqueness proof; it does not show that a solver can find the support or that no distant
% alternative exists. Geometries as diag_geo, diag_geo2, diag_geo3: G1 Section VI full great circles with both poles (16 rotations /
% 180 deg), G2 half arc pole to pole (32 / 360), G3 full ring with 15-deg polar gaps (16 / 180), G4 half arc 5-175 (32 / 360), G5 one-sided
% arch 15-165 as built (32 / 360). Writes <results>/jac_rot_geometry.mat and prints the table. One thread, about 0.5 GB.
%   T = jac_rot_geometry()            seed 1, tol = [0.999 0.9999]
if nargin < 1 || isempty(seed), seed = 1; end
if nargin < 2 || isempty(tols), tols = [0.999 0.9999]; end
cd(fileparts(fileparts(fileparts(mfilename('fullpath'))))); setup; cpr_settings('paper'); maxNumCompThreads(1); warning('off', 'all');
vFile = fullfile(cpr_paths().data, 'fullwave', 'fwc_validation.mat'); Sv = load(vFile); V = Sv.V;
modes = V.A.modes; nf = 16; L = 12; nm = size(modes, 1); a = V.A.a(:, 1:nf); f = V.freq(1:nf).';
fac = (1i .^ (modes(:, 1) + (modes(:, 3) == 1))).';                          % far-field convention (diag_geo)
geos = {};
ts = (0:31).' * 11.25; far = ts > 180; th = ts; th(far) = 360 - ts(far); ph = zeros(size(ts)); ph(far) = pi;
sg = ones(size(ts)); sg(far) = -1; th = min(max(th, 0.01), 179.99) * pi/180;
[Yt, Yp] = vsh_eval(L, th, ph); geos(end+1, :) = {'G1 full circles with both poles', [sg .* Yt; Yp] .* fac, (0:15) * 11.25}; %#ok<AGROW>
th2 = linspace(0.01, 179.99, 16).' * pi/180; [Yt, Yp] = vsh_eval(L, th2, zeros(16, 1));
geos(end+1, :) = {'G2 half arc pole to pole', [Yt; Yp] .* fac, (0:31) * 11.25};
th3 = [(15:10:165) (15:10:165)].' * pi/180; ph3 = [zeros(1, 16) pi*ones(1, 16)].'; [Yt, Yp] = vsh_eval(L, th3, ph3);
geos(end+1, :) = {'G3 full ring with 15-deg polar gaps', [Yt; Yp] .* fac, (0:15) * 11.25};
th4 = linspace(5, 175, 16).' * pi/180; [Yt, Yp] = vsh_eval(L, th4, zeros(16, 1));
geos(end+1, :) = {'G4 half arc 5-175 deg', [Yt; Yp] .* fac, (0:31) * 11.25};
th5 = (15:10:165).' * pi/180; [Yt, Yp] = vsh_eval(L, th5, zeros(16, 1));
geos(end+1, :) = {'G5 one-sided arch 15-165 deg', [Yt; Yp] .* fac, (0:31) * 11.25};
T = struct('geometry', {}, 'tol', {}, 'nSupport', {}, 'm', {}, 'JsNull', {}, 'rotMarginDb', {}, 'rotCountLt30', {}, 'rotDirs', {}, ...
    'rotFreqMarginDb', {}, 'rotFreqCountLt30', {}, 'rotFreqDirs', {});
for gi = 1:size(geos, 1)
    R = geos{gi, 2}; psi = geos{gi, 3}; np = size(R, 1); nrot = numel(psi);
    E = exp(-1i * modes(:, 2).' .* (psi(:) * pi/180));                        % nrot x nm
    G = cell(1, nf); Xm = zeros(np, nrot, nf);
    for k = 1:nf
        Gk = zeros(np*nrot, nm);
        for l = 1:nrot, Gk((l-1)*np + (1:np), :) = R .* E(l, :); end
        G{k} = Gk; Xm(:, :, k) = reshape(Gk * a(:, k), np, nrot);
    end
    c = cpr_cfg('arch', 'IC', 'snrdB', 100, 'nd', nf, 'K', 4, 'ratio', 0.3, 'vshL', L); c.nlegs = 0;
    Dm = struct('Xring', Xm, 'f', f(:).', 'nlegsAvail', 0, 'Xleg', [], 'anchor', 1, 'legShift', 1);
    S = make_readings_fw(c, Dm, seed, 'ideal');
    x = Xm(:); z = S.Phi * x; u = z ./ abs(z);
    [ii, jj, vv] = find(S.Phi);
    rot = floor(mod(jj - 1, np*nrot) / np) + 1; frq = floor((jj - 1) / (np*nrot)) + 1;
    dz = 1i * vv .* x(jj);                                                    % d z_i / d g for the entry's rotation
    cu = real(conj(u(ii)) .* dz);
    JG1 = full(sparse(ii, rot, cu, S.m, nrot));
    JG2 = full(sparse(ii, (frq - 1)*nrot + rot, cu, S.m, nrot*nf));
    s = a(:); [~, o] = sort(abs(s).^2, 'descend');
    for tol = tols
        kk = find(cumsum(abs(s(o)).^2) >= tol * norm(s)^2, 1); sup = o(1:kk); [qm, fk] = ind2sub([nm nf], sup);
        rows = {}; cols = {}; vals = {};
        for t = 1:kk
            col = G{fk(t)}(:, qm(t)); rr = (fk(t) - 1)*np*nrot + (1:np*nrot).';
            rows{t} = rr; cols{t} = t*ones(size(rr)); vals{t} = col; %#ok<AGROW>
        end
        Psi = sparse(vertcat(rows{:}), vertcat(cols{:}), vertcat(vals{:}), np*nrot*nf, kk);
        H = S.Phi * Psi; uH = full(conj(u) .* H); Js = [real(uH), -imag(uH)];
        [Us, Ss, ~] = svd(Js, 'econ'); ss = diag(Ss); nul = sum(ss < ss(1)*1e-9); Us = Us(:, ss > ss(1)*1e-10);
        [m1, c1, d1] = margin(JG1, Us, ones(1, nrot));
        [m2, c2, d2] = margin(JG2, Us, ones(1, nrot*nf));
        T(end+1) = struct('geometry', geos{gi, 1}, 'tol', tol, 'nSupport', kk, 'm', S.m, 'JsNull', nul, 'rotMarginDb', m1, ...
            'rotCountLt30', c1, 'rotDirs', d1, 'rotFreqMarginDb', m2, 'rotFreqCountLt30', c2, 'rotFreqDirs', d2); %#ok<AGROW>
        fprintf('%-38s tol %.2f%% (%4d coeffs, m = %d): Js null dim %d | per rotation: smallest %7.1f dB, %d of %d below -30 dB | per rotation and frequency: smallest %7.1f dB, %d of %d\n', ...
            geos{gi, 1}, 100*tol, kk, S.m, nul, m1, c1, d1, m2, c2, d2);
    end
end
save(fullfile(cpr_paths().results, 'jac_rot_geometry.mat'), 'T', 'seed');
end

function [mdb, cnt, nd] = margin(JG, Us, gauge)
% smallest normalized singular value (dB) of (I - Us Us') JG Q, Q spanning the complement of the gauge (global phase)
Q = null(gauge); W = JG * Q; Wo = W - Us * (Us' * W); [~, Rw] = qr(W, 0);
sv = svd(Wo / Rw); sv = max(sv, realmin); d = sort(20*log10(sv)); mdb = d(1); cnt = sum(d < -30); nd = numel(d);
end
