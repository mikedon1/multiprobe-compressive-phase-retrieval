function T = fwc_ring_true(V, srcDir)
%FWC_RING_TRUE  Measured AUT-to-probe transmissions of the TRUE full ring (65-port HFSS model, 16 rotations over 180 deg).
%   T = fwc_ring_true(V, srcDir)
% V       fwc_validate_hfss result (its map V.map lists the 32 arch ports as (pol, ring_index) in the order of V.R)
% srcDir  folder with sweep_fullring_rot<k>_psi<psi>.s65p (one per rotation) and fullring_port_map.csv
% T.Xring 64 x nrot x nf: rows 1..32 the near half in the order of V.map, rows 33..64 the far half in the same order, the layout
%         of the emulated ring of fwc_run_trials ('Ring', 'full'), so the same recovery operator rows apply
% T.psi   rotation angles (rad), T.freq, T.files, T.mapFull (the full-ring port map)
% The far half is the near half turned by 180 deg about the rotation axis, so a far port is paired with the near port of the same
% (polarization, ring index) and the far responses of the model follow from the near ones by the factor (-1)^m (fwc_run_trials).
d = dir(fullfile(srcDir, 'sweep_fullring_rot*.s65p')); [~, i] = sort({d.name}); d = d(i);
if isempty(d), error('fwc_ring_true: no sweep_fullring_rot*.s65p files in %s', srcDir); end
Mf = readtable(fullfile(srcDir, 'fullring_port_map.csv'), 'TextType', 'string');
np = height(V.map);
jn = zeros(np, 1); jf = zeros(np, 1);
for k = 1:np
    a = find(Mf.half == "near" & Mf.pol == V.map.pol(k) & Mf.ring_index == V.map.ring_index(k));
    b = find(Mf.half == "far" & Mf.pol == V.map.pol(k) & Mf.ring_index == V.map.ring_index(k));
    if numel(a) ~= 1 || numel(b) ~= 1, error('fwc_ring_true: port %s %d not unique in the full-ring map', V.map.pol(k), V.map.ring_index(k)); end
    jn(k) = Mf.touchstone_index(a); jf(k) = Mf.touchstone_index(b);
end
nrot = numel(d); T.psi = zeros(1, nrot); T.files = {d.name};
for r = 1:nrot
    [S, f, names, ~, vars] = read_touchstone(fullfile(d(r).folder, d(r).name));
    if size(S, 1) ~= 65, error('fwc_ring_true: %s has %d ports, expected 65', d(r).name, size(S, 1)); end
    if r == 1, T.freq = f(:).'; T.Xring = zeros(2*np, nrot, numel(f)); T.names = names; end
    if any(abs(f(:).' - T.freq) > 1) || ~isequal(names, T.names), error('fwc_ring_true: %s differs from the first file', d(r).name); end
    T.psi(r) = sscanf(vars.AUT_psi, '%f') * pi/180;
    T.Xring(1:np, r, :) = permute(S(jn, 1, :), [1 3 2]);                  % AUT (port 1) -> near-half probe ports
    T.Xring(np+1:2*np, r, :) = permute(S(jf, 1, :), [1 3 2]);             % AUT -> far-half probe ports
end
if any(abs(T.freq - V.freq(:).') > 1), error('fwc_ring_true: frequencies differ from the validation data'); end
T.mapFull = Mf;
end
