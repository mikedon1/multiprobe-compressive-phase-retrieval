function [R, info] = fwc_probe_resp_all(TsAll, radii, Lreg, Lout, rq)
%FWC_PROBE_RESP_ALL  fwc_probe_resp for every port at once (round 7, 2026-10-20): the same mathematics, one factorization
%per frequency instead of one per port and frequency.
%   [R, info] = fwc_probe_resp_all(TsAll, radii, Lreg, Lout, rq)
% TsAll{p}{s}: read_aedt_table near-field export of port p on sphere s (radii(s)); all ports share the sampling grid
% and the frequency list (checked). For every frequency the regular-wave basis B at the sample points of both spheres
% does not depend on the port, so b = B \ Y solves all ports as the columns of Y with ONE factorization (mldivide
% makes the same rank decisions as for a single column, since they depend on B only). The quadrature-sphere bases are
% port independent as well.
% R (np x nmode x nf): R(p, q, f) as fwc_probe_resp(...) returns for port p, R(q, f).
% info{p}: .fitres, .radres (nf x 1), .modes, .freq, as fwc_probe_resp.
c0 = 299792458;
np = numel(TsAll); ns = numel(radii);
freq = TsAll{1}{1}.freq(:).'; nf = numel(freq);
[xg, wg] = gauss_leg(2*Lreg + 4); thq = acos(xg); nphq = 2*Lreg + 8; phq = (0:nphq-1) * 2*pi/nphq;
[TQ, PQ] = ndgrid(thq, phq); WQ = wg(:) * ones(1, nphq) * (2*pi/nphq);
wq = WQ(:) * rq^2;
nmo = 2*Lout*(Lout+2);
% sample grid per sphere from port 1; every port must share it
samp = cell(1, ns);
for s = 1:ns
    T = TsAll{1}{s};
    if any(abs(T.freq(:).' - freq) > 1), error('fwc_probe_resp_all:freq', 'frequency lists differ'); end
    [TH, PH] = ndgrid(T.x(:) * pi/180, T.phi(:).' * pi/180);
    th = TH(:); ph = PH(:);
    n = [sin(th).*cos(ph), sin(th).*sin(ph), cos(th)];
    [~, keep] = unique(round(n * 1e8), 'rows', 'stable');
    samp{s}.th = th(keep); samp{s}.ph = ph(keep); samp{s}.keep = keep;
    samp{s}.w = sqrt(max(sin(th(keep)), 1e-3));
    samp{s}.E = cell(np, 3);
    for p = 1:np
        Tp = TsAll{p}{s};
        if ~isequal(Tp.x(:), T.x(:)) || ~isequal(Tp.phi(:), T.phi(:)) || any(abs(Tp.freq(:).' - freq) > 1)
            error('fwc_probe_resp_all:grid', 'port %d sphere %d has a different grid or frequency list', p, s);
        end
        samp{s}.E{p, 1} = reshape(Tp.c.NearEX, [], nf); samp{s}.E{p, 2} = reshape(Tp.c.NearEY, [], nf);
        samp{s}.E{p, 3} = reshape(Tp.c.NearEZ, [], nf);
    end
end
R = zeros(np, nmo, nf);
fitres = zeros(np, nf); radres = nan(np, nf);
for f = 1:nf
    k = 2*pi*freq(f) / c0;
    Bs = cell(1, ns); Ys = cell(1, ns);
    for s = 1:ns
        S = samp{s};
        Et = zeros(numel(S.th), np); Ep = Et;
        for p = 1:np
            E = [S.E{p, 1}(S.keep, f), S.E{p, 2}(S.keep, f), S.E{p, 3}(S.keep, f)];
            [tp, pp, Et(:, p), Ep(:, p)] = fwc_hfss2rot(S.th, S.ph, E);
        end
        [Rt, Rp] = sph_modes(Lreg, k, radii(s), tp, pp, 'reg');     % directions are port independent
        Bs{s} = [S.w .* Rt; S.w .* Rp]; Ys{s} = [S.w .* Et; S.w .* Ep];
    end
    B = vertcat(Bs{:}); Y = vertcat(Ys{:});
    b = B \ Y;                                                      % all ports, one factorization
    fitres(:, f) = (vecnorm(B*b - Y) ./ vecnorm(Y)).';
    if ns == 2
        if rcond(Bs{1}'*Bs{1}) > 1e-14 && rcond(Bs{2}'*Bs{2}) > 1e-14
            b1 = Bs{1} \ Ys{1}; b2 = Bs{2} \ Ys{2};
            radres(:, f) = (vecnorm(b1 - b2) ./ vecnorm(b)).';
        end
    end
    clear B Y Bs Ys
    [Ebt, Ebp, Hbt, Hbp] = sph_modes(Lreg, k, rq, TQ(:), PQ(:), 'reg');
    ebtA = Ebt*b; ebpA = Ebp*b; hbtA = Hbt*b; hbpA = Hbp*b;
    [Ot, Op, Oht, Ohp, modes] = sph_modes(Lout, k, rq, TQ(:), PQ(:), 'out');
    for p = 1:np
        R(p, :, f) = (((Ot .* hbpA(:, p) - Op .* hbtA(:, p)) - (ebtA(:, p) .* Ohp - ebpA(:, p) .* Oht)).' * wq).';
    end
    fprintf('fwc_probe_resp_all: frequency %d of %d done\n', f, nf);
end
info = cell(1, np);
for p = 1:np
    info{p} = struct('fitres', fitres(p, :).', 'radres', radres(p, :).', 'modes', modes, 'freq', freq);
end
end

function [x, w] = gauss_leg(n)
b = 0.5 ./ sqrt(1 - (2*(1:n-1)).^(-2));
[V, D] = eig(diag(b, 1) + diag(b, -1));
x = diag(D); w = 2 * V(1, :).'.^2;
end
