function fwc_pack_ring(hfssDir, nfDirs, outDir)
%FWC_PACK_RING  Pack the raw HFSS exports of the TRUE full ring (65-port model) into compressed MATLAB files for release.
%   fwc_pack_ring(hfssDir, nfDirs, outDir)
% hfssDir  folder with sweep_fullring_rot*.s65p (16 rotations over 180 deg), fullring_port_map.csv, ringalone_fullring_map.csv
% nfDirs   folder or cell array of folders with ringalone_fullring_nf/NF12_<terminal>.csv and NF17_<terminal>.csv (ring alone, near fields on
%          spheres of 12 and 17 cm, one port excited at a time); a port found in several folders is taken from the first
% outDir   default <repository>/data/hfss_fullring, receives
%   hfss_ring_rotations.mat   S (65 x 65 x freq x 16) single precision, psi (deg), port names, variables, the full-ring port map
%   hfss_ring_nf_r12.mat,     E (theta x phi x freq x 3 x port) in V/m, the ports with both sphere files in the order of ringalone_fullring_map.csv
%   hfss_ring_nf_r17.mat      (map stored with them)
% Complex values are single precision (about -144 dB, far below the solver accuracy). fwc_ring_true and fwc_check_ring_resp read these files
% when the raw exports are absent.
if nargin < 3 || isempty(outDir), outDir = fullfile(cpr_paths().data, 'hfss_fullring'); end
if ischar(nfDirs) || isstring(nfDirs), nfDirs = {char(nfDirs)}; end
if ~exist(outDir, 'dir'), mkdir(outDir); end
sgl = @(z) single(z);
d = dir(fullfile(hfssDir, 'sweep_fullring_rot*.s65p')); [~, i] = sort({d.name}); d = d(i);
rot.psi = zeros(1, numel(d)); rot.vars = cell(1, numel(d));
for k = 1:numel(d)
    [S, f, names, ~, vars] = read_touchstone(fullfile(d(k).folder, d(k).name));
    if k == 1, rot.S = complex(zeros([size(S) numel(d)], 'single')); rot.freq = f(:).'; rot.portNames = names; end
    if any(abs(f(:).' - rot.freq) > 1) || ~isequal(names, rot.portNames), error('fwc_pack_ring: %s differs', d(k).name); end
    rot.S(:, :, :, k) = sgl(S); rot.psi(k) = sscanf(vars.AUT_psi, '%f'); rot.vars{k} = vars;
end
rot.files = {d.name};
rot.mapFull = readtable(fullfile(hfssDir, 'fullring_port_map.csv'), 'TextType', 'string');
rot.units = 'S: [port x port x freq x rotation], Touchstone port order (1 = AUT), 50 ohm; psi in deg; freq in Hz';
save(fullfile(outDir, 'hfss_ring_rotations.mat'), '-struct', 'rot', '-v7');
M = readtable(fullfile(hfssDir, 'ringalone_fullring_map.csv'), 'TextType', 'string');
names = {'NF12', 'NF17'}; rad = [0.12 0.17]; src = strings(height(M), 2);
for p = 1:height(M), for s = 1:2
    for q = 1:numel(nfDirs)
        f = fullfile(nfDirs{q}, 'ringalone_fullring_nf', sprintf('%s_%s.csv', names{s}, M.nf_terminal(p)));
        if isfile(f) && dir(f).bytes > 2^20, src(p, s) = f; break, end
    end
end, end
have = all(src ~= "", 2); M = M(have, :); src = src(have, :);
fprintf('ports with both sphere files: %d of 64\n', height(M));
for s = 1:2
    nf = struct('radius', rad(s), 'terminal', M.nf_terminal.', 'map', M);
    for p = 1:height(M)
        T = read_aedt_table(char(src(p, s)));
        if p == 1
            nf.theta = T.x(:).'; nf.phi = T.phi(:).'; nf.freq = T.freq(:).';
            nf.E = complex(zeros([size(T.c.NearEX) 3 height(M)], 'single'));
        elseif ~isequal(T.x(:).', nf.theta) || ~isequal(T.phi(:).', nf.phi) || any(abs(T.freq(:).' - nf.freq) > 1)
            error('fwc_pack_ring: %s has a different grid', src(p, s));
        end
        nf.E(:, :, :, 1, p) = sgl(T.c.NearEX); nf.E(:, :, :, 2, p) = sgl(T.c.NearEY); nf.E(:, :, :, 3, p) = sgl(T.c.NearEZ);
    end
    nf.units = 'E: [theta x phi x freq x (Ex, Ey, Ez) x port] in V/m, ports in the order of map; deg; Hz';
    save(fullfile(outDir, sprintf('hfss_ring_nf_r%d.mat', round(rad(s)*100))), '-struct', 'nf', '-v7');
end
dd = dir(fullfile(outDir, 'hfss_ring_*.mat'));
for k = 1:numel(dd), fprintf('%-24s %7.1f MB\n', dd(k).name, dd(k).bytes/2^20); end
end
