function T = fwc_run_trials(seeds, varargin)
%FWC_RUN_TRIALS  Section VII trials on the full-arc HFSS data with embedded-probe calibration (2026-10-19).
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
p.parse(varargin{:});
over = struct2nv(p.Unmatched);
V = p.Results.Validation; if isempty(V), V = fwc_validate_hfss(); end
D = struct('Xring', V.Smeas, 'R', V.R, 'modes', V.A.modes, 'psi', V.psi(:).' * pi/180, 'f', V.freq(:).', ...
    'aTrue', V.A.a);
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
c = cpr_cfg('arch', 'IC', 'snrdB', 50, over{:});
fprintf('Section VII trials: %d ports x %d rotations x %d frequencies, m/n %.2f, SNR %d dB, P %d, K %d\n', ...
    size(D.Xring, 1), size(D.Xring, 2), size(D.Xring, 3), c.ratio, c.snrdB, c.P, c.K);
T = struct('seed', {}, 'nmse', {}, 'time', {}, 'numTry', {}, 'aHat', {}, 'perfreqFactor', {});
for s = seeds
    r = fwc_recover(c, s, D);
    T(end+1) = struct('seed', s, 'nmse', r.nmse, 'time', r.time, 'numTry', r.info.numTry, 'aHat', r.aHat, ...
        'perfreqFactor', r.perfreqFactor); %#ok<AGROW>
    fprintf('  seed %2d: per-frequency %.1f dB, constant+delay %.1f dB, global phase %.1f dB, %d tries, %.0f s\n', ...
        s, r.nmse.perfreqScale, r.nmse.globalDelay, r.nmse.globalPhase, r.info.numTry, r.time);
end
v = arrayfun(@(t) t.nmse.perfreqScale, T); g = arrayfun(@(t) t.nmse.globalDelay, T);
fprintf('per-frequency score: mean %.1f dB, best %.1f, worst %.1f; constant+delay: mean %.1f, best %.1f, worst %.1f\n', ...
    mean(v), min(v), max(v), mean(g), min(g), max(g));
tag = ''; if o.PerturbAzDeg > 0 || isfinite(o.PerturbRelDb), tag = sprintf('_az%g_rel%g', o.PerturbAzDeg, o.PerturbRelDb); end
if strcmp(o.ProbeModel, 'free'), tag = ['_free' tag]; end
save(fullfile(cpr_paths().results, ['fwc_trials' tag '.mat']), 'T', 'c', 'o');
end

function nv = struct2nv(s)
f = fieldnames(s); nv = cell(1, 2*numel(f));
for i = 1:numel(f), nv{2*i-1} = f{i}; nv{2*i} = s.(f{i}); end
end
