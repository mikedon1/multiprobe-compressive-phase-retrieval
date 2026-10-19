function V = fwc_validate_hfss(cfg)
%FWC_VALIDATE_HFSS  Step 4 of sim/hfss/DUALPOL_RING_PLAN.md: predict the solved dual-ring rotations
%from the AUT-alone far field (Solve A) and the ring-alone probe near fields (Solve B).
%   V = fwc_validate_hfss(cfg)
% cfg fields (defaults for the private-repo layout):
%   (defaults: the FULL-ARC model of 2026-10-19 -- offset AUT, ring 15-165 deg; the old 20-90 ring used
%   autalone_farfield.csv, ringalone_dual_*, sweep_dualpol_rot*)
%   ffFile   Solve A far-field CSV
%   nfDir    Solve B near-field CSVs                  (results/ringalone_dual_nf)
%   mapFile  Solve B port map CSV: columns nf_terminal, full_port (Touchstone index in the full model,
%            from geometry), written by the Solve B script
%   sGlob    full-model rotation files                (results/sweep_dualpol_rot*.s33p)
%   radii    sphere radii (m) and names               ([0.12 0.17], {'NF12','NF17'})
%   L        AUT mode degree (12); Lreg regular-wave degree (20); rq quadrature radius (0.145 m)
% V: nmse_db (2 x nf, rows = rotation sense +1 / -1), kappa, perPort (sense +1), A (Solve A fit),
%    R (probe responses), info per port, and the measured/predicted arrays.
R0 = 'C:\Users\miked\Projects\cpr-multiprobe-antenna\sim\hfss\results';
def = struct('ffFile', fullfile(R0, 'autalone_offset_farfield.csv'), 'nfDir', fullfile(R0, 'ringalone_fullarc_nf'), ...
    'mapFile', fullfile(R0, 'ringalone_fullarc_map.csv'), 'sGlob', fullfile(R0, 'sweep_fullarc_rot*.s33p'), ...
    'radii', [0.12 0.17], 'nfNames', {{'NF12', 'NF17'}}, 'L', 12, 'Lreg', 20, 'rq', 0.145);
if nargin < 1, cfg = struct(); end
for f = fieldnames(def).', if ~isfield(cfg, f{1}), cfg.(f{1}) = def.(f{1}); end, end

A = fwc_aut_modes(read_aedt_table(cfg.ffFile), cfg.L);
fprintf('Solve A: fit residual %.1f to %.1f dB; energy above degree 10: %.1e (worst frequency)\n', ...
    20*log10(min(A.relres)), 20*log10(max(A.relres)), max(sum(A.degE(:, 11:end), 2)));

d = dir(cfg.sGlob); files = fullfile({d.folder}, {d.name});
[Sall, freq, ~, ~, ~] = read_touchstone(files{1}); %#ok<ASGLU>
nf = numel(freq);
if any(abs(freq(:).' - A.freq) > 1), error('fwc_validate_hfss: Solve A and sweep frequencies differ'); end
psi = zeros(1, numel(files)); Sm = [];
for k = 1:numel(files)
    [S, ~, ~, ~, vars] = read_touchstone(files{k});
    psi(k) = sscanf(vars.AUT_psi, '%f');
    Sm(:, k, :) = S(2:end, 1, :); %#ok<AGROW>   % AUT (port 1) -> probe ports 2..33
end

M = readtable(cfg.mapFile, 'TextType', 'string');
have = true(height(M), 1);                          % use only ports whose sphere files all exist (partial runs)
for p = 1:height(M)
    for s = 1:numel(cfg.radii)
        have(p) = have(p) && isfile(fullfile(cfg.nfDir, sprintf('%s_%s.csv', cfg.nfNames{s}, M.nf_terminal(p))));
    end
end
M = M(have, :);
fprintf('ports with near-field data: %d of %d\n', height(M), numel(have));
np = height(M); R = zeros(np, size(A.a, 1), nf); info = cell(1, np);
for p = 1:np
    Ts = cell(1, numel(cfg.radii));
    for s = 1:numel(cfg.radii)
        Ts{s} = read_aedt_table(fullfile(cfg.nfDir, sprintf('%s_%s.csv', cfg.nfNames{s}, M.nf_terminal(p))));
    end
    [Rp, info{p}] = fwc_probe_resp(Ts, cfg.radii, cfg.Lreg, cfg.L, cfg.rq);
    R(p, :, :) = permute(Rp, [3 1 2]);
end
fprintf('Solve B: regular-wave fit residual worst %.1f dB; radius consistency worst %.1f dB\n', ...
    20*log10(max(cellfun(@(x) max(x.fitres), info))), 20*log10(max(cellfun(@(x) max(x.radres), info))));
Smeas = Sm(M.full_port - 1, :, :);              % rows in map order
V.nmse_db = zeros(2, nf); V.kappa = zeros(2, nf);
senses = [1 -1];
for k = 1:2
    Sp = fwc_predict(A.a, A.modes, R, psi*pi/180, senses(k));
    [V.nmse_db(k, :), V.kappa(k, :), pp] = fwc_fit_scale(Sp, Smeas);
    if k == 1, V.perPort = pp; V.Spred = Sp; end
    fprintf('sense %+d: NMSE %.1f .. %.1f dB over frequencies\n', senses(k), min(V.nmse_db(k, :)), max(V.nmse_db(k, :)));
end
V.Smeas = Smeas; V.psi = psi; V.freq = A.freq; V.A = A; V.R = R; V.info = info; V.map = M;
end
