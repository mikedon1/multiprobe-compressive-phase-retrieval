function out = fwc_check_nf_port(file, r, pos, ax)
%FWC_CHECK_NF_PORT  Check that a near-field export belongs to the probe that was meant to be driven.
%   out = fwc_check_nf_port(file, r, pos, ax)
% file: AEDT near-field data table (re/im NearEX/EY/EZ, X = Theta, columns per Phi and Freq, HFSS frame)
% on a sphere of radius r (m) around the origin; pos: probe centre (m, HFSS frame); ax: probe wire axis.
% For each frequency:
%   polFrac  fraction of the field energy along the probe axis (near 1 for a dipole at ~1 m: its field
%            near the origin is polarized along the wire)
%   dirErr   angle (deg) between the probe direction pos/|pos| and the beam-scan maximum
%            u* = argmax_u |sum_n (E_n . ax) exp(-1i k u . r_n)| (the field arrives from the probe)
c0 = 299792458;
T = read_aedt_table(file);
[TH, PH] = ndgrid(T.x(:) * pi/180, T.phi(:).' * pi/180);
n = [sin(TH(:)).*cos(PH(:)), sin(TH(:)).*sin(PH(:)), cos(TH(:))];
[~, keep] = unique(round(n * 1e8), 'rows', 'stable'); n = n(keep, :); P = r * n;
w = max(sin(acos(n(:, 3))), 1e-3);                     % area weight of the equal-angle grid
ax = ax(:).' / norm(ax); u0 = pos(:).' / norm(pos);
[tu, pu] = ndgrid((0:1:180) * pi/180, (0:1:359) * pi/180);
U = [sin(tu(:)).*cos(pu(:)), sin(tu(:)).*sin(pu(:)), cos(tu(:))];
nf = numel(T.freq); out.polFrac = zeros(1, nf); out.dirErr = zeros(1, nf); out.freq = T.freq;
for f = 1:nf
    E = [reshape(T.c.NearEX(:, :, f), [], 1), reshape(T.c.NearEY(:, :, f), [], 1), reshape(T.c.NearEZ(:, :, f), [], 1)];
    E = E(keep, :);
    e = E * ax.';
    out.polFrac(f) = sum(w .* abs(e).^2) / sum(w .* sum(abs(E).^2, 2));
    k = 2*pi*T.freq(f) / c0;
    S = abs(exp(-1i * k * (U * P.')) * (w .* e));
    [~, i] = max(S);
    out.dirErr(f) = acosd(max(-1, min(1, U(i, :) * u0.')));
    out.uMax(f, :) = U(i, :);
end
end
