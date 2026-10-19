function [R, info] = fwc_probe_resp(Ts, radii, Lreg, Lout, rq)
%FWC_PROBE_RESP  Embedded receiving response of one probe port to every outgoing AUT mode (Solve B).
%   [R, info] = fwc_probe_resp(Ts, radii, Lreg, Lout, rq)
% Ts: cell of read_aedt_table near-field exports of this port excited (ring alone, AUT absent), one
% per sphere radius, each with X = Theta (deg), columns per Phi and Freq, quantities re/im of NearEX,
% NearEY, NearEZ (HFSS frame). radii: sphere radii (m), same order. Lreg: degree of the regular-wave
% (standing-wave) fit of the probe field, well above the AUT degree (about k*r_max + 8). Lout: AUT mode
% degree. rq: radius (m) of the quadrature sphere of the reciprocity integral (between the AUT's
% minimum sphere and the probes; the result does not depend on it if the fit is good).
% R (nmode(Lout) x nf): R(q, f) = integral over the sphere of (E_q x H_p - E_p x H_q) . r-hat, with
% E_q, H_q the outgoing mode q (sph_modes 'out') and E_p, H_p the probe field from its regular-wave fit
% (sph_modes 'reg'; H from the fitted coefficients, so no H export is needed). Then the port response
% to an AUT with mode coefficients a is sum_q a_q R_q, up to one constant per frequency
% (validated with exact dipole fields in sph_modes_test.m).
% info.fitres (nf x 1): relative residual of the regular-wave fit on the two spheres;
% info.radres (nf x 1): relative mismatch between the coefficients fitted on each sphere alone (if
% both alone are well posed; NaN otherwise) -- an internal consistency check.
c0 = 299792458;
freq = Ts{1}.freq(:).'; nf = numel(freq);
[xg, wg] = gauss_leg(2*Lreg + 4); thq = acos(xg); nphq = 2*Lreg + 8; phq = (0:nphq-1) * 2*pi/nphq;
[TQ, PQ] = ndgrid(thq, phq); WQ = wg(:) * ones(1, nphq) * (2*pi/nphq);
wq = WQ(:) * rq^2;
R = zeros(2*Lout*(Lout+2), nf); info.fitres = zeros(nf, 1); info.radres = nan(nf, 1);
samp = cell(1, numel(Ts));
for s = 1:numel(Ts)
    T = Ts{s};
    if any(abs(T.freq(:).' - freq) > 1), error('fwc_probe_resp:freq', 'frequency lists differ'); end
    [TH, PH] = ndgrid(T.x(:) * pi/180, T.phi(:).' * pi/180);
    th = TH(:); ph = PH(:);
    n = [sin(th).*cos(ph), sin(th).*sin(ph), cos(th)];
    [~, keep] = unique(round(n * 1e8), 'rows', 'stable');
    samp{s}.th = th(keep); samp{s}.ph = ph(keep); samp{s}.keep = keep;
    samp{s}.w = sqrt(max(sin(th(keep)), 1e-3));
    samp{s}.E = {reshape(T.c.NearEX, [], nf), reshape(T.c.NearEY, [], nf), reshape(T.c.NearEZ, [], nf)};
end
for f = 1:nf
    k = 2*pi*freq(f) / c0;
    B = []; y = []; Bs = cell(1, numel(Ts)); ys = Bs;
    for s = 1:numel(Ts)
        S = samp{s};
        E = [S.E{1}(S.keep, f), S.E{2}(S.keep, f), S.E{3}(S.keep, f)];
        [tp, pp, Et, Ep] = fwc_hfss2rot(S.th, S.ph, E);
        [Rt, Rp] = sph_modes(Lreg, k, radii(s), tp, pp, 'reg');
        Bs{s} = [S.w .* Rt; S.w .* Rp]; ys{s} = [S.w .* Et; S.w .* Ep];
        B = [B; Bs{s}]; y = [y; ys{s}]; %#ok<AGROW>
    end
    b = B \ y;
    info.fitres(f) = norm(B*b - y) / norm(y);
    if numel(Ts) == 2
        b1 = Bs{1} \ ys{1}; b2 = Bs{2} \ ys{2};
        if rcond(Bs{1}'*Bs{1}) > 1e-14 && rcond(Bs{2}'*Bs{2}) > 1e-14
            info.radres(f) = norm(b1 - b2) / norm(b);
        end
    end
    [Ebt, Ebp, Hbt, Hbp] = sph_modes(Lreg, k, rq, TQ(:), PQ(:), 'reg');
    ebt = Ebt*b; ebp = Ebp*b; hbt = Hbt*b; hbp = Hbp*b;
    [Ot, Op, Oht, Ohp, modes] = sph_modes(Lout, k, rq, TQ(:), PQ(:), 'out');
    R(:, f) = ((Ot .* hbp - Op .* hbt) - (ebt .* Ohp - ebp .* Oht)).' * wq;
end
info.modes = modes; info.freq = freq;
end

function [x, w] = gauss_leg(n)
b = 0.5 ./ sqrt(1 - (2*(1:n-1)).^(-2));
[V, D] = eig(diag(b, 1) + diag(b, -1));
x = diag(D); w = 2 * V(1, :).'.^2;
end
