function V = fwc_validate_hfss(cfg)
%FWC_VALIDATE_HFSS  Step 4 of sim/hfss/DUALPOL_RING_PLAN.md: predict the solved dual-ring rotations
%from the AUT-alone far field (Solve A) and the ring-alone probe near fields (Solve B).
%   V = fwc_validate_hfss(cfg)
% cfg fields (defaults: the raw HFSS exports in <repository>/data/hfss_fullarc, or in the folder named by the
% environment variable CPR_HFSS):
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
R0 = getenv('CPR_HFSS'); if isempty(R0), R0 = fullfile(cpr_paths().data, 'hfss_fullarc'); end
def = struct('ffFile', fullfile(R0, 'autalone_offset_farfield.csv'), 'nfDir', fullfile(R0, 'ringalone_fullarc_nf'), ...
    'mapFile', fullfile(R0, 'ringalone_fullarc_map.csv'), 'sGlob', fullfile(R0, 'sweep_fullarc_rot*.s33p'), ...
    'radii', [0.12 0.17], 'nfNames', {{'NF12', 'NF17'}}, 'L', 12, 'Lreg', 20, 'rq', 0.145);
if nargin < 1, cfg = struct(); end
for f = fieldnames(def).', if ~isfield(cfg, f{1}), cfg.(f{1}) = def.(f{1}); end, end
% packed release files (fwc_pack_hfss): used when the raw exports are absent, or with cfg.packed = true
if ~isfield(cfg, 'packed'), cfg.packed = ~isfile(cfg.ffFile) && isfile(fullfile(R0, 'hfss_aut.mat')); end
if cfg.packed, V = validate_packed(R0, cfg); return, end

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
if ~isfield(cfg, 'fast'), cfg.fast = true; end            % round 7: one factorization per frequency for all ports
if cfg.fast
    TsAll = cell(1, np);
    for p = 1:np
        TsAll{p} = cell(1, numel(cfg.radii));
        for s = 1:numel(cfg.radii)
            TsAll{p}{s} = read_aedt_table(fullfile(cfg.nfDir, sprintf('%s_%s.csv', cfg.nfNames{s}, M.nf_terminal(p))));
        end
        if mod(p, 8) == 0, fprintf('read the near-field files of %d of %d ports\n', p, np); end
    end
    [R, info] = fwc_probe_resp_all(TsAll, cfg.radii, cfg.Lreg, cfg.L, cfg.rq);
    clear TsAll
else
    for p = 1:np
        Ts = cell(1, numel(cfg.radii));
        for s = 1:numel(cfg.radii)
            Ts{s} = read_aedt_table(fullfile(cfg.nfDir, sprintf('%s_%s.csv', cfg.nfNames{s}, M.nf_terminal(p))));
        end
        [Rp, info{p}] = fwc_probe_resp(Ts, cfg.radii, cfg.Lreg, cfg.L, cfg.rq);
        R(p, :, :) = permute(Rp, [3 1 2]);
    end
end
V = finish(A, R, info, Sm, psi, M);
end

function V = validate_packed(R0, cfg)
% the same computation from the packed files of fwc_pack_hfss (single precision converted to double)
Pa = load(fullfile(R0, 'hfss_aut.mat'));
T = struct('x', Pa.theta(:), 'phi', Pa.phi(:).', 'freq', Pa.freq(:).');
T.c.rETheta = double(Pa.rETheta); T.c.rEPhi = double(Pa.rEPhi);
A = fwc_aut_modes(T, cfg.L);
fprintf('Solve A (packed): fit residual %.1f to %.1f dB\n', 20*log10(min(A.relres)), 20*log10(max(A.relres)));
Pr = load(fullfile(R0, 'hfss_rotations.mat'));
if any(abs(Pr.freq - A.freq) > 1), error('fwc_validate_hfss: Solve A and sweep frequencies differ'); end
Sm = permute(double(Pr.S(2:end, 1, :, :)), [1 4 3 2]);   % AUT (port 1) -> probe ports 2..33: [port x rot x freq]
psi = Pr.psi; M = Pr.mapRing;
files = {'hfss_nf_r12.mat', 'hfss_nf_r17.mat'};
np = height(M); TsAll = cell(1, np);
for s = 1:numel(files)
    Pn = load(fullfile(R0, files{s}));
    if abs(Pn.radius - cfg.radii(s)) > 1e-9 || ~isequal(cellstr(Pn.terminal), cellstr(M.nf_terminal.'))
        error('fwc_validate_hfss: %s does not match the radii or the port map', files{s});
    end
    for p = 1:np
        Tp = struct('x', Pn.theta(:), 'phi', Pn.phi(:).', 'freq', Pn.freq(:).');
        Tp.c.NearEX = double(Pn.E(:, :, :, 1, p)); Tp.c.NearEY = double(Pn.E(:, :, :, 2, p));
        Tp.c.NearEZ = double(Pn.E(:, :, :, 3, p));
        TsAll{p}{s} = Tp;
    end
end
[R, info] = fwc_probe_resp_all(TsAll, cfg.radii, cfg.Lreg, cfg.L, cfg.rq);
V = finish(A, R, info, Sm, psi, M);
end

function V = finish(A, R, info, Sm, psi, M)
nf = numel(A.freq);
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
