function [Yt, Yp, modes] = vsh_eval(L, th, ph)
%VSH_EVAL  Tangential vector spherical harmonics (far-field spherical-wave modes) up to degree L.
%   [Yt, Yp, modes] = vsh_eval(L, th, ph)
% th, ph: directions (rad, column vectors). Yt, Yp (npts x ncols): theta and phi components of
% each mode. Mode (l, m, type), l = 1..L, m = -l..l:
%   type 1 (M):  r x grad Y_lm = -(i m / sin th) Y_lm theta-hat + dY_lm/dth phi-hat
%   type 2 (N):  grad Y_lm     =  dY_lm/dth theta-hat + (i m / sin th) Y_lm phi-hat
% each divided by sqrt(l(l+1)), so the modes are orthonormal over the sphere. Y_lm uses MATLAB's
% normalized associated Legendre functions (Condon-Shortley phase). modes = [l m type] per column.
% A pole direction (sin th < 1e-6) is evaluated at th = 1e-4 (or pi - 1e-4) along its ph: the limit
% along the meridian, as a probe on a scan circle through the pole sees it.
th = th(:); ph = ph(:);
pole = sin(th) < 1e-6;
th(pole & th < 1) = 1e-4; th(pole & th >= 1) = pi - 1e-4;
d = 1e-6; st = sin(th);
nm = 2*L*(L+2); Yt = zeros(numel(th), nm); Yp = Yt; modes = zeros(nm, 3); q = 0;
for l = 1:L
    P0 = legendre(l, cos(th).', 'norm').';                 % npts x (l+1), m = 0..l
    dP = (legendre(l, cos(th + d).', 'norm').' - legendre(l, cos(th - d).', 'norm').') / (2*d);
    g = 1/sqrt(l*(l+1)) / sqrt(2*pi);
    for m = -l:l
        sg = 1; if m < 0, sg = (-1)^m; end
        e = exp(1i*m*ph);
        p0 = sg*P0(:, abs(m)+1).*e; dp = sg*dP(:, abs(m)+1).*e;
        q = q + 1; Yt(:, q) = g*(-1i*m*p0./st); Yp(:, q) = g*dp;     modes(q, :) = [l m 1];
        q = q + 1; Yt(:, q) = g*dp;             Yp(:, q) = g*(1i*m*p0./st); modes(q, :) = [l m 2];
    end
end
end
