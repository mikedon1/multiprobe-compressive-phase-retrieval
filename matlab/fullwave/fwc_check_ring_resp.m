function out = fwc_check_ring_resp(V, srcDir)
%FWC_CHECK_RING_RESP  Embedded probe responses of the FULL ring (near-field exports of step87) against the arch-only calibration.
%   out = fwc_check_ring_resp(V, srcDir)
% V       fwc_validate_hfss result (V.R: 32 arch ports x modes x frequencies; V.map: their (pol, ring_index))
% srcDir  folder with ringalone_fullring_map.csv and ringalone_fullring_nf/NF12_<terminal>.csv, NF17_<terminal>.csv
% For every port with both sphere files, the response R_p(l, m, f) is computed with the same code as the arch (fwc_probe_resp_all, regular-wave
% degree 20, output degree 12, quadrature radius 14.5 cm) and compared with
%   near-half ports: the arch response of the same (pol, ring index), which is the calibration used so far;
%   far-half ports : (-1)^m times the arch response of the same (pol, ring index) (the 180-deg symmetry used for the far half).
% Reported per frequency (dB): the relative difference, and the difference after one complex gain per port and frequency (shape only).
% A difference well below the model error of Section VII (about -31 to -38 dB) means that the arch-only calibration is adequate for the
% full ring and the full 64-port Solve B is not needed.
radii = [0.12 0.17]; names = {'NF12', 'NF17'}; Lreg = 20; L = 12; rq = 0.145;
packed = isfile(fullfile(srcDir, 'hfss_ring_nf_r12.mat')) && ~isfolder(fullfile(srcDir, 'ringalone_fullring_nf'));   % files of fwc_pack_ring
if packed
    Pn = {load(fullfile(srcDir, 'hfss_ring_nf_r12.mat')), load(fullfile(srcDir, 'hfss_ring_nf_r17.mat'))};
    if abs(Pn{1}.radius - radii(1)) > 1e-9 || abs(Pn{2}.radius - radii(2)) > 1e-9, error('fwc_check_ring_resp: packed radii differ'); end
    M = Pn{1}.map; np = height(M); TsAll = cell(1, np);
    for p = 1:np
        for s = 1:2
            Tp = struct('x', Pn{s}.theta(:), 'phi', Pn{s}.phi(:).', 'freq', Pn{s}.freq(:).');
            Tp.c.NearEX = double(Pn{s}.E(:, :, :, 1, p)); Tp.c.NearEY = double(Pn{s}.E(:, :, :, 2, p)); Tp.c.NearEZ = double(Pn{s}.E(:, :, :, 3, p));
            TsAll{p}{s} = Tp;
        end
    end
    fprintf('ports with near-field data (packed): %d\n', np);
else
    M = readtable(fullfile(srcDir, 'ringalone_fullring_map.csv'), 'TextType', 'string');
    nfDir = fullfile(srcDir, 'ringalone_fullring_nf');
    have = false(height(M), 1);
    for p = 1:height(M), have(p) = all(arrayfun(@(s) isfile(fullfile(nfDir, sprintf('%s_%s.csv', names{s}, M.nf_terminal(p)))), 1:2)); end
    M = M(have, :); np = height(M);
    if np == 0, error('fwc_check_ring_resp: no port with both sphere files in %s', nfDir); end
    fprintf('ports with near-field data: %d\n', np);
    TsAll = cell(1, np);
    for p = 1:np
        TsAll{p} = cell(1, 2);
        for s = 1:2, TsAll{p}{s} = read_aedt_table(fullfile(nfDir, sprintf('%s_%s.csv', names{s}, M.nf_terminal(p)))); end
    end
end
[Rn, info] = fwc_probe_resp_all(TsAll, radii, Lreg, L, rq);       % np x nmode x nf
fprintf('regular-wave fit residual worst %.1f dB\n', 20*log10(max(cellfun(@(x) max(x.fitres), info))));
modes = V.A.modes; sgn = reshape((-1) .^ modes(:, 2), 1, []);
nf = size(Rn, 3); out.table = cell(np, 1); out.diff = zeros(np, nf); out.shape = zeros(np, nf);
fprintf('%-10s %-5s %-5s %3s | %-26s | %-26s\n', 'terminal', 'half', 'pol', 'idx', 'difference (dB): min..max', 'shape only (dB): min..max');
for p = 1:np
    k = find(V.map.pol == M.pol(p) & V.map.ring_index == M.ring_index(p));
    if numel(k) ~= 1, error('fwc_check_ring_resp: arch port %s %d not unique', M.pol(p), M.ring_index(p)); end
    Rref = squeeze(V.R(k, :, :));                                       % nmode x nf
    if M.half(p) == "far", Rref = Rref .* sgn.'; end
    Rp = squeeze(Rn(p, :, :));
    for q = 1:nf
        a = Rref(:, q); b = Rp(:, q);
        out.diff(p, q) = 20*log10(norm(b - a) / norm(a));
        g = (a' * b) / (a' * a); out.shape(p, q) = 20*log10(norm(b - g*a) / norm(b));
    end
    fprintf('%-10s %-5s %-5s %3d | %8.1f ..%8.1f        | %8.1f ..%8.1f\n', M.nf_terminal(p), M.half(p), M.pol(p), M.ring_index(p), ...
        min(out.diff(p, :)), max(out.diff(p, :)), min(out.shape(p, :)), max(out.shape(p, :)));
end
out.map = M; out.R = Rn;
fprintf('Section VII model error for comparison: -38.3 .. -30.8 dB per frequency\n');
end
