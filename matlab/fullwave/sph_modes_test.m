function sph_modes_test()
%SPH_MODES_TEST  Validate sph_modes + regular-wave fit + reciprocity coupling with Hertzian dipoles (exact fields).
maxNumCompThreads(1);
addpath(genpath('C:\Users\miked\Projects\multiprobe-compressive-phase-retrieval\matlab'));
k = 2*pi/0.1;                                     % 3 GHz
Lout = 14; Lreg = 22;
% quadrature sphere for the reciprocity integral
[xg, wg] = gauss_leg(48); thq = acos(xg); phq = (0:95)*2*pi/96;
[TQ, PQ] = ndgrid(thq, phq); WQ = wg(:) * ones(1, 96) * (2*pi/96);
% sample spheres for the probe-field fit (5 deg grid, as the HFSS export)
[TS, PS] = ndgrid((0:5:180)*pi/180, (0:5:355)*pi/180);
rng(3);
res = [];
for trial = 1:6
    ra = (rand(3,1) - 0.5) * 0.08;               % AUT dipole within 4 cm of the origin
    pa = randn(3,1) + 1i*randn(3,1);
    rb = randn(3,1); rb = rb / norm(rb) * 1.0;    % probe dipole at 1 m
    pb = randn(3,1);
    % AUT mode coefficients from its far field (as from the HFSS AUT-alone export)
    [TF, PF] = ndgrid((1:2:179)*pi/180, (0:4:356)*pi/180);
    n = sphdir(TF(:), PF(:));
    F = cross(cross(n, repmat(pa.', size(n,1), 1), 2), n, 2) .* exp(1i*k*(n*ra));
    [Ft, Fp] = cart2tp(F, TF(:), PF(:));
    [Yt, Yp, modes] = vsh_eval(Lout, TF(:), PF(:));
    w = sqrt(sin(TF(:)));
    c = [w.*Yt; w.*Yp] \ [w.*Ft; w.*Fp];
    l = modes(:, 1); t = modes(:, 3);
    a = c ./ (1i.^(l + (t == 1)));               % far field: TE ~ j^(l+1), TM ~ j^l
    % probe field near the origin: exact, sampled on two spheres, fitted with regular modes
    B = []; y = [];
    for rs = [0.12 0.17]
        P = rs * sphdir(TS(:), PS(:));
        [E, ~] = dipole_field(pb, rb, P, k);
        [Et, Ep] = cart2tp(E, TS(:), PS(:));
        [Rt, Rp] = sph_modes(Lreg, k, rs, TS(:), PS(:), 'reg');
        B = [B; Rt; Rp]; y = [y; Et; Ep]; %#ok<AGROW>
    end
    b = B \ y; fitErr = norm(B*b - y)/norm(y);
    % H check on the quadrature sphere: fitted vs exact
    rq = 0.15; P = rq * sphdir(TQ(:), PQ(:));
    [Ex, Hx] = dipole_field(pb, rb, P, k);
    [Ebt, Ebp] = cart2tp(Ex, TQ(:), PQ(:)); [Hbt, Hbp] = cart2tp(Hx, TQ(:), PQ(:));
    [~, ~, Rht, Rhp] = sph_modes(Lreg, k, rq, TQ(:), PQ(:), 'reg');
    hErr = norm([Rht*b - Hbt; Rhp*b - Hbp]) / norm([Hbt; Hbp]);
    % reciprocity coupling of each outgoing mode with the probe field (fitted b)
    [Rt2, Rp2, Rh2t, Rh2p] = sph_modes(Lreg, k, rq, TQ(:), PQ(:), 'reg');
    ebt = Rt2*b; ebp = Rp2*b; hbt = Rh2t*b; hbp = Rh2p*b;
    [Ot, Op, Oht, Ohp] = sph_modes(Lout, k, rq, TQ(:), PQ(:), 'out');
    wq = WQ(:) * rq^2;
    Rcoef = ((Ot .* hbp - Op .* hbt) - (ebt .* Ohp - ebp .* Oht)).' * wq;   % one per AUT mode
    pred = sum(a .* Rcoef);
    % exact coupling: p_b . E_a(r_b)
    [Ea, ~] = dipole_field(pa, ra, rb.', k);
    truth = Ea * pb;
    res(end+1, :) = [fitErr, hErr, pred/truth]; %#ok<AGROW>
    fprintf('trial %d: regular fit %.1e, H from fit vs exact %.1e, pred/truth = %.6f %+.6fi\n', ...
        trial, fitErr, hErr, real(pred/truth), imag(pred/truth));
end
q = res(:, 3); fprintf('spread of pred/truth: %.2e (relative)\n', max(abs(q - mean(q)))/abs(mean(q)));
end

function n = sphdir(th, ph)
n = [sin(th).*cos(ph), sin(th).*sin(ph), cos(th)];
end

function [Ft, Fp] = cart2tp(F, th, ph)
tt = [cos(th).*cos(ph), cos(th).*sin(ph), -sin(th)];
pp = [-sin(ph), cos(ph), zeros(size(ph))];
Ft = sum(F .* tt, 2); Fp = sum(F .* pp, 2);
end

function [E, H] = dipole_field(p, r0, P, k)
% Hertzian dipole p at r0, fields at rows of P (eps = mu = c = 1, time e^{jwt})
R = P - r0(:).'; d = sqrt(sum(R.^2, 2)); n = R ./ d;
g = exp(-1i*k*d);
pn = n * p(:);
nxp = cross(n, repmat(p(:).', size(n,1), 1), 2);
E = (1/(4*pi)) * (k^2 * cross(nxp, n, 2) .* (g./d) + (3*n.*pn - p(:).') .* ((1./d.^3 + 1i*k./d.^2) .* g));
H = (k^2/(4*pi)) * nxp .* (g./d) .* (1 - 1i./(k*d));
end

function [x, w] = gauss_leg(n)
b = 0.5 ./ sqrt(1 - (2*(1:n-1)).^(-2));
[V, D] = eig(diag(b, 1) + diag(b, -1));
x = diag(D); w = 2 * V(1, :).'.^2;
end
