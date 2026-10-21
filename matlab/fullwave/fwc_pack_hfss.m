function fwc_pack_hfss(srcDir, outDir)
%FWC_PACK_HFSS  Pack the raw HFSS exports of the Section VII model into compressed MATLAB files for release.
%   fwc_pack_hfss(srcDir, outDir)
% srcDir holds the AEDT exports (autalone_offset_farfield.csv, autalone_offset_rot00.s1p, sweep_fullarc_rot*.s33p,
% ringalone_fullarc_rot00.s32p, ringalone_fullarc_map.csv, fullarc_port_map.csv, ringalone_fullarc_nf/NF12_*.csv and
% NF17_*.csv). outDir (default <repository>/data/hfss_fullarc) receives:
%   hfss_aut.mat         Solve A, AUT alone: far field rE_theta, rE_phi (V) on a 2 x 5 deg grid, and its S11
%   hfss_rotations.mat   the full 33-port model at 32 rotations: S (33 x 33 x 16 x 32), psi, port names, variables;
%                        the ring alone at psi = 0: Sring (32 x 32 x 16); both port maps
%   hfss_nf_r12.mat,     Solve B, ring alone: E_x, E_y, E_z (V/m) on spheres of radius 12 and 17 cm, 5 x 5 deg grid,
%   hfss_nf_r17.mat      for each of the 32 probe ports excited in turn (the other ports matched)
% Complex values are stored in single precision (relative rounding ~6e-8, about -144 dB, far below the solver's
% accuracy of about 0.02 dB); the exports carry 15 digits. fwc_validate_hfss reads these files when the CSV and
% Touchstone exports are absent (cfg.packed).
if nargin < 2 || isempty(outDir), outDir = fullfile(cpr_paths().data, 'hfss_fullarc'); end
if ~exist(outDir, 'dir'), mkdir(outDir); end
sgl = @(z) single(z);
note = ['HFSS-IE 2016.2 exports of the Section VII model (dual-polarized arch 15-165 deg, 16 positions, 33 ports, ' ...
    'offset tilted dipole AUT), packed by fwc_pack_hfss. Complex values in single precision.'];

% ---- Solve A ----
T = read_aedt_table(fullfile(srcDir, 'autalone_offset_farfield.csv'));
aut.theta = T.x(:).'; aut.phi = T.phi(:).'; aut.freq = T.freq(:).';
aut.rETheta = sgl(T.c.rETheta); aut.rEPhi = sgl(T.c.rEPhi);       % [theta x phi x freq], V
aut.units = 'rETheta, rEPhi in V (r times E), theta/phi in deg, freq in Hz';
[S1, f1] = read_touchstone(fullfile(srcDir, 'autalone_offset_rot00.s1p'));
aut.S11 = sgl(squeeze(S1)); aut.S11freq = f1(:).';
aut.note = note;
save(fullfile(outDir, 'hfss_aut.mat'), '-struct', 'aut', '-v7');

% ---- full model, 32 rotations, and the ring alone ----
d = dir(fullfile(srcDir, 'sweep_fullarc_rot*.s33p'));
[~, i] = sort({d.name}); d = d(i);
rot.psi = zeros(1, numel(d)); rot.vars = cell(1, numel(d));
for k = 1:numel(d)
    [S, f, names, ~, vars] = read_touchstone(fullfile(d(k).folder, d(k).name));
    if k == 1, rot.S = complex(zeros([size(S) numel(d)], 'single')); rot.freq = f(:).'; rot.portNames = names; end
    if any(abs(f(:).' - rot.freq) > 1) || ~isequal(names, rot.portNames), error('fwc_pack_hfss: %s differs', d(k).name); end
    rot.S(:, :, :, k) = sgl(S);
    rot.psi(k) = sscanf(vars.AUT_psi, '%f'); rot.vars{k} = vars;
end
rot.files = {d.name};
[Sr, fr, nr] = read_touchstone(fullfile(srcDir, 'ringalone_fullarc_rot00.s32p'));
rot.Sring = sgl(Sr); rot.SringFreq = fr(:).'; rot.SringPortNames = nr;
rot.mapFull = readtable(fullfile(srcDir, 'fullarc_port_map.csv'), 'TextType', 'string');
rot.mapRing = readtable(fullfile(srcDir, 'ringalone_fullarc_map.csv'), 'TextType', 'string');
rot.units = 'S: [port x port x freq x rotation], Touchstone port order (1 = AUT), 50 ohm; psi in deg; freq in Hz';
rot.note = note;
save(fullfile(outDir, 'hfss_rotations.mat'), '-struct', 'rot', '-v7');

% ---- Solve B near fields ----
M = rot.mapRing;
for s = {{'NF12', 0.12, 'hfss_nf_r12.mat'}, {'NF17', 0.17, 'hfss_nf_r17.mat'}}
    nm = s{1}{1}; nf = struct('radius', s{1}{2}, 'terminal', M.nf_terminal.');
    for p = 1:height(M)
        T = read_aedt_table(fullfile(srcDir, 'ringalone_fullarc_nf', sprintf('%s_%s.csv', nm, M.nf_terminal(p))));
        if p == 1
            nf.theta = T.x(:).'; nf.phi = T.phi(:).'; nf.freq = T.freq(:).';
            nf.E = complex(zeros([size(T.c.NearEX) 3 height(M)], 'single'));
        elseif ~isequal(T.x(:).', nf.theta) || ~isequal(T.phi(:).', nf.phi) || any(abs(T.freq(:).' - nf.freq) > 1)
            error('fwc_pack_hfss: %s port %s has a different grid', nm, M.nf_terminal(p));
        end
        nf.E(:, :, :, 1, p) = sgl(T.c.NearEX); nf.E(:, :, :, 2, p) = sgl(T.c.NearEY); nf.E(:, :, :, 3, p) = sgl(T.c.NearEZ);
        if mod(p, 8) == 0, fprintf('%s: %d of %d ports\n', nm, p, height(M)); end
    end
    nf.units = 'E: [theta x phi x freq x (Ex, Ey, Ez) x port] in V/m, ports in the order of terminal (map order); deg; Hz';
    nf.note = note;
    save(fullfile(outDir, s{1}{3}), '-struct', 'nf', '-v7');
end
dd = dir(fullfile(outDir, 'hfss_*.mat'));
for k = 1:numel(dd), fprintf('%-22s %7.1f MB\n', dd(k).name, dd(k).bytes/2^20); end
end
