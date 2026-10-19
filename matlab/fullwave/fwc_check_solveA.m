% OBSOLETE (2026-10-19): written for the OLD 20-90 deg ring (v7: ports 2..17 / 33..18, theta = 20..90). It hard-codes
% that port order and arc and would compute the wrong thing on the full-arc model (fullarc_port_map.csv). Kept for the
% record of the 2026-10-18 check; use fwc_validate_hfss for the full-arc model.
function C = fwc_check_solveA(ffFile, sGlob, L)
error('fwc_check_solveA:obsolete', 'written for the old 20-90 deg ring; see the header');
%FWC_CHECK_SOLVEA  Check Solve A and the rotation sense on the solved dual-ring rotations, without
%Solve B: the free-space AUT far field (VSH fit), rotated by psi, predicts each probe's transmission
%up to one complex gain per probe and frequency (an ideal-probe model, T5's "one gain per probe").
%   C = fwc_check_solveA(ffFile, sGlob, L)
% Ring (sim/hfss/DUALPOL_RING_PLAN.md): rotation frame z' = +y, x' = +z, y' = +x; probes on phi' = 0 at
% theta' = 20 + i*70/15; ports 2..17 theta-hat (co = E_theta'), ports 18..33 phi-hat with theta
% DESCENDING (cross = E_phi'). Rotation psi multiplies mode (l, m) by exp(-1i*sense*m*psi).
% Reports, per sense and polarization, the NMSE (dB) with one gain per probe and frequency, and with one
% gain per frequency only (no probe calibration at all).
R0 = 'C:\Users\miked\Projects\cpr-multiprobe-antenna\sim\hfss\results';
if nargin < 1 || isempty(ffFile), ffFile = fullfile(R0, 'autalone_farfield.csv'); end
if nargin < 2 || isempty(sGlob), sGlob = fullfile(R0, 'sweep_dualpol_rot*.s33p'); end
if nargin < 3, L = 12; end
T = read_aedt_table(ffFile);
fprintf('far-field table: %d theta x %d phi x %d freq, unit scales %s\n', numel(T.x), numel(T.phi), numel(T.freq), mat2str(T.unitScales));
A = fwc_aut_modes(T, L);
A16 = fwc_aut_modes(T, 16);
fprintf('Solve A VSH fit residual, L = %d: %.1f .. %.1f dB; L = 16: %.1f .. %.1f dB\n', L, ...
    20*log10(min(A.relres)), 20*log10(max(A.relres)), 20*log10(min(A16.relres)), 20*log10(max(A16.relres)));
fprintf('energy in degrees 11-12 (L = 12 fit), worst frequency: %.1e; degrees 13-16 (L = 16 fit): %.1e\n', ...
    max(sum(A.degE(:, 11:12), 2)), max(sum(A16.degE(:, 13:16), 2)));
d = dir(sGlob); files = fullfile({d.folder}, {d.name});
psi = zeros(1, numel(files)); Sm = [];
for k = 1:numel(files)
    [S, f, ~, ~, vars] = read_touchstone(files{k});
    psi(k) = sscanf(vars.AUT_psi, '%f');
    Sm(:, k, :) = S(2:33, 1, :); %#ok<AGROW>
end
if any(abs(f(:).' - A.freq) > 1), error('frequency lists differ'); end
thp = (20 + (0:15) * 70/15) * pi/180;
[Yt, Yp, modes] = vsh_eval(L, thp(:), zeros(16, 1));
co = Sm(1:16, :, :); cr = Sm(32:-1:17, :, :);          % cross rows reordered to theta ascending
m = modes(:, 2); nf = numel(A.freq);
C.psi = psi; C.relres = A.relres; C.degE = A.degE;
for sense = [1 -1]
    for pol = 1:2
        if pol == 1, Y = Yt; meas = co; nm = 'co (E_theta'')'; else, Y = Yp; meas = cr; nm = 'cross (E_phi'')'; end
        pred = zeros(size(meas));
        for fi = 1:nf
            pred(:, :, fi) = Y * (A.c(:, fi) .* exp(-1i * sense * m * psi(:).' * pi/180));
        end
        % one complex gain per probe and frequency
        e1 = 0; e0 = 0; e2 = 0;
        for fi = 1:nf
            P = pred(:, :, fi); M = meas(:, :, fi);
            g = sum(conj(P) .* M, 2) ./ sum(abs(P).^2, 2);
            e1 = e1 + sum(abs(g .* P - M).^2, 'all'); e0 = e0 + sum(abs(M).^2, 'all');
            k = (P(:)' * M(:)) / (P(:)' * P(:));
            e2 = e2 + sum(abs(k * P - M).^2, 'all');
        end
        fprintf('sense %+d, %-14s: gain per probe+freq %.1f dB, gain per freq only %.1f dB  (%d rotations)\n', ...
            sense, nm, 10*log10(e1/e0), 10*log10(e2/e0), numel(psi));
        C.nmse(sense == [1 -1], pol, :) = [10*log10(e1/e0), 10*log10(e2/e0)];
    end
end
end
