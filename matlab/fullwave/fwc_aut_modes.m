function A = fwc_aut_modes(T, L)
%FWC_AUT_MODES  Free-space AUT mode coefficients from the HFSS AUT-alone far field (Solve A).
%   A = fwc_aut_modes(T, L)
% T: read_aedt_table of the far-field export (X = Theta in deg, columns per Phi and Freq, quantities
% re/im of rETheta and rEPhi, HFSS frame). L: truncation degree of the fit.
% The field is moved to the rotation-axis frame (fwc_hfss2rot), duplicate directions (theta < 0 with
% phi, theta = +/-180) are removed, and the VSH coefficients c are fitted per frequency by weighted least
% squares (weight sqrt(sin theta), the area element of the equal-angle grid). The mode coefficients of
% the outgoing waves of sph_modes.m follow as a = c ./ j^(l+1) (TE) and c ./ j^l (TM), up to one
% constant per frequency (k and the e^{-jkr}/r normalization), which the calibration absorbs.
% A.a (nmode x nf), A.c, A.modes, A.freq, A.relres (fit residual per frequency), A.degE (energy per
% degree, nf x L), A.tp, A.pp, A.Et, A.Ep (the samples used, primed frame).
[TH, PH] = ndgrid(T.x(:) * pi/180, T.phi(:).' * pi/180);
nf = numel(T.freq);
Ft = reshape(T.c.rETheta, numel(TH), nf); Fp = reshape(T.c.rEPhi, numel(TH), nf);
th = TH(:); ph = PH(:);
% HFSS spherical -> Cartesian, then unique directions
n = [sin(th).*cos(ph), sin(th).*sin(ph), cos(th)];
[~, keep] = unique(round(n * 1e8), 'rows', 'stable');
th = th(keep); ph = ph(keep); Ft = Ft(keep, :); Fp = Fp(keep, :);
tt = [cos(th).*cos(ph), cos(th).*sin(ph), -sin(th)];
pv = [-sin(ph), cos(ph), zeros(size(ph))];
[Yt, Yp, modes] = deal([]);
A.c = zeros(2*L*(L+2), nf); A.relres = zeros(1, nf);
w = sqrt(max(abs(sin(th)), 1e-3));                % HFSS-grid area weight (poles kept once, small weight)
for f = 1:nf
    E = Ft(:, f) .* tt + Fp(:, f) .* pv;
    [tp, pp, Et, Ep] = fwc_hfss2rot(th, ph, E);
    if f == 1
        [Yt, Yp, modes] = vsh_eval(L, tp, pp);
        B = [w .* Yt; w .* Yp];
        [Q, R] = qr(B, 0);
    end
    y = [w .* Et; w .* Ep];
    A.c(:, f) = R \ (Q' * y);
    A.relres(f) = norm(B * A.c(:, f) - y) / norm(y);
    A.Et(:, f) = Et; A.Ep(:, f) = Ep;
end
l = modes(:, 1); t = modes(:, 3);
A.a = A.c ./ (1i .^ (l + (t == 1)));
A.modes = modes; A.freq = T.freq(:).'; A.tp = tp; A.pp = pp; A.L = L;
A.degE = zeros(nf, L);
for d = 1:L
    A.degE(:, d) = sum(abs(A.c(l == d, :)).^2, 1).';
end
A.degE = A.degE ./ sum(A.degE, 2);
end
