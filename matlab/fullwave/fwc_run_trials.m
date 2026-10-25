function T = fwc_run_trials(seeds, varargin)
%FWC_RUN_TRIALS  Section VII trials on the HFSS data with embedded-probe calibration (2026-10-19; emulated full ring by
%default since 2026-10-21, option 'Ring').
%   T = fwc_run_trials(1:10)                         IC, m/n 0.3, 50 dB SNR (the paper's Section VII)
%   T = fwc_run_trials(1, 'maxTry', 1, 'nitEM', 2)   extra name/value pairs override cpr_cfg
%   T = fwc_run_trials(1:10, 'Validation', V)        reuse a fwc_validate_hfss result instead of recomputing it
% Data: fwc_validate_hfss() (Solve A modes, Solve B embedded responses of every port with near-field data, the
% measured AUT->probe responses of the 32 rotations, ports in map order). Each trial: fwc_recover (on/off
% switches, P ports, K tones, random cables, noise, PR-GAMP with the calibrated operator). Scores per trial:
% perfreqScale (main), perfreqPhase, globalPhase, globalScale, globalDelay (one constant + one delay), model.
% Saves <results>/fwc_trials.mat (T: one entry per seed) and prints a summary (mean of dB, best, worst).
p = inputParser; p.KeepUnmatched = true; p.addParameter('Validation', []);
p.addParameter('PerturbAzDeg', 0);    % robustness (2026-10-19): rms azimuthal position error of each probe position in
                                      %   the CALIBRATION (both ports of a position share it), deg; the data stay HFSS
p.addParameter('PerturbRelDb', -inf); % random complex error of every embedded-response entry, rms relative to the rms
                                      %   response of its port, dB (e.g. -40: the probe-equalization residual of [Varela 2024])
p.addParameter('PerturbSeed', 1);
p.addParameter('ProbeModel', 'embedded'); % 'embedded' (the paper) | 'free': ideal point probes, no coupling (fwc_free_resp)
p.addParameter('OutTag', '');         % appended to the output file name (parallel workers on disjoint seeds; merge after)
p.addParameter('RingDir', fullfile(cpr_paths().data, 'hfss_fullring'));   % 'true' ring: folder with sweep_fullring_rot*.s65p + fullring_port_map.csv
p.addParameter('Ring', 'full');       % 'true' (2026-10-22): the TRUE full ring from the 65-port HFSS model (see RingDir);
                                      % 'full' (the paper, 2026-10-21): the arch at psi + 180 deg is the far half of a full
                                      %   ring at psi, so the 32 computed rotations form a 64-port ring at 16 rotations
                                      %   (coupling between the two halves is not in the single-arch HFSS model);
                                      %   'arc': the one-sided arch at 32 rotations (fails for the dipole AUT, see
                                      %   results_vsh/drivers/diag_geo.m and the supplement's sampling-geometry table)
p.addParameter('RingData', []);       % 'true' ring: struct with Xring, psi (rad) [, Rring] instead of the HFSS files (mom_make_ring)
p.addParameter('RingCal', 'arch');    % 'true' ring calibration: 'arch' (arch response + 180-deg symmetry) | 'ring' (RingData.Rring, ring alone)
p.addParameter('Polish', false);      % true (2026-10-23): refine each PR-GAMP estimate with fwc_polish (held-out stopping), default off
p.parse(varargin{:});
over = struct2nv(p.Unmatched);
V = p.Results.Validation;
vFile = fullfile(cpr_paths().data, 'fullwave', 'fwc_validation.mat');   % processed HFSS data released with the code
if isempty(V) && exist(vFile, 'file'), S = load(vFile); V = S.V; end
if isempty(V), V = fwc_validate_hfss(); end                              % from the raw exports (data/hfss_fullarc)
D = struct('Xring', V.Smeas, 'R', V.R, 'modes', V.A.modes, 'psi', V.psi(:).' * pi/180, 'f', V.freq(:).', ...
    'aTrue', V.A.a, 'polish', p.Results.Polish);
o = p.Results;
if strcmp(o.ProbeModel, 'free')
    D.R = fwc_free_resp(V);
    fprintf('probe model: FREE SPACE (ideal point probes, no coupling)\n');
end
if o.PerturbAzDeg > 0 || isfinite(o.PerturbRelDb)
    rs = RandStream('mt19937ar', 'Seed', 7000 + o.PerturbSeed);
    [np, nm, nf] = size(D.R);
    if o.PerturbAzDeg > 0                     % a probe moved by delta in phi = the AUT turned by -delta: e^{+i m delta}
        ri = V.map.ring_index(:); [u, ~, iu] = unique(ri);
        dPos = (o.PerturbAzDeg*pi/180) * randn(rs, numel(u), 1); dP = dPos(iu);
        D.R = D.R .* exp(1i * dP .* reshape(D.modes(:, 2), 1, nm));
    end
    if isfinite(o.PerturbRelDb)
        sc = sqrt(mean(abs(D.R).^2, [2 3]));  % rms response per port
        E = (randn(rs, np, nm, nf) + 1i*randn(rs, np, nm, nf)) / sqrt(2);
        D.R = D.R + 10^(o.PerturbRelDb/20) * sc .* E;
    end
    fprintf('calibration perturbed: azimuth %.2f deg rms, entries %.0f dB\n', o.PerturbAzDeg, o.PerturbRelDb);
end
if strcmp(o.Ring, 'full')                     % emulated full ring: ports [arch; arch turned by 180 deg], 16 rotations
    nr = numel(D.psi); h = nr/2;
    if mod(nr, 2) || any(abs(D.psi(h+1:end) - D.psi(1:h) - pi) > 1e-9)
        error('fwc_run_trials: Ring full needs rotations paired 180 deg apart');
    end
    D.Xring = [D.Xring(:, 1:h, :); D.Xring(:, h+1:end, :)];
    D.R = [D.R; D.R .* reshape((-1) .^ D.modes(:, 2), 1, [])];   % e^{-i m (psi + pi)} = (-1)^m e^{-i m psi}
    D.psi = D.psi(1:h);
elseif strcmp(o.Ring, 'true')                 % TRUE full ring (65-port HFSS model): measured far-half ports, 16 rotations over 180 deg
    if ~isempty(o.RingData), Tr = o.RingData; else, Tr = fwc_ring_true(V, o.RingDir); end   % RingData: e.g. mom_make_ring (moment-method ring)
    if numel(Tr.psi) ~= 16 || any(abs(Tr.psi - (0:15) * pi/16) > 1e-6), error('fwc_run_trials: Ring true needs the 16 rotations psi = 0:11.25:168.75 deg'); end
    D.Xring = Tr.Xring; D.psi = Tr.psi;
    D.R = [D.R; D.R .* reshape((-1) .^ D.modes(:, 2), 1, [])];   % far-half responses from the 180-deg symmetry (calibration: arch only)
    if strcmp(o.RingCal, 'ring'), D.R = Tr.Rring; end              % calibration of the full ring alone (RingData with field Rring)
elseif ~strcmp(o.Ring, 'arc')
    error('fwc_run_trials: unknown Ring %s', o.Ring);
end
c = cpr_cfg('arch', 'IC', 'snrdB', 50, over{:});
fprintf('Section VII trials: %d ports x %d rotations x %d frequencies, m/n %.2f, SNR %d dB, P %d, K %d\n', ...
    size(D.Xring, 1), size(D.Xring, 2), size(D.Xring, 3), c.ratio, c.snrdB, c.P, c.K);
T = struct('seed', {}, 'nmse', {}, 'time', {}, 'numTry', {}, 'aHat', {}, 'perfreqFactor', {}, 'nres_dB', {}, 'log', {});
for s = seeds
    r = fwc_recover(c, s, D);
    T(end+1) = struct('seed', s, 'nmse', r.nmse, 'time', r.time, 'numTry', r.info.numTry, 'aHat', r.aHat, ...
        'perfreqFactor', r.perfreqFactor, 'nres_dB', r.info.nres_dB, 'log', r.info.log); %#ok<AGROW>   % accepted residual and per-restart log (review item 6.2)
    fprintf('  seed %2d: per-frequency %.1f dB, constant+delay %.1f dB, global phase %.1f dB, %d tries, %.0f s\n', ...
        s, r.nmse.perfreqScale, r.nmse.globalDelay, r.nmse.globalPhase, r.info.numTry, r.time);
end
v = arrayfun(@(t) t.nmse.perfreqScale, T); g = arrayfun(@(t) t.nmse.globalDelay, T);
fprintf('per-frequency score: mean %.1f dB, best %.1f, worst %.1f; constant+delay: mean %.1f, best %.1f, worst %.1f\n', ...
    mean(v), min(v), max(v), mean(g), min(g), max(g));
tag = ''; if o.PerturbAzDeg > 0 || isfinite(o.PerturbRelDb), tag = sprintf('_az%g_rel%g', o.PerturbAzDeg, o.PerturbRelDb); end
if strcmp(o.ProbeModel, 'free'), tag = ['_free' tag]; end
if strcmp(o.Ring, 'arc'), tag = ['_arc' tag]; end
if strcmp(o.Ring, 'true'), tag = ['_true' tag]; end
tag = [tag o.OutTag];
save(fullfile(cpr_paths().results, ['fwc_trials' tag '.mat']), 'T', 'c', 'o');
end

function nv = struct2nv(s)
f = fieldnames(s); nv = cell(1, 2*numel(f));
for i = 1:numel(f), nv{2*i-1} = f{i}; nv{2*i} = s.(f{i}); end
end
