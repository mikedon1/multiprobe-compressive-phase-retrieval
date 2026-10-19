function [Et, Ep, Ht, Hp, modes] = sph_modes(L, k, r, th, ph, kind, eta)
%SPH_MODES  Tangential E and H of the spherical-wave modes up to degree L (time convention e^{jwt}).
%   [Et, Ep, Ht, Hp, modes] = sph_modes(L, k, r, th, ph, kind, eta)
% r, th, ph: points (vectors of equal length or r scalar). kind: 'out' (spherical Hankel h_l^(2),
% outgoing) or 'reg' (spherical Bessel j_l, standing; fields regular at the origin). eta: wave
% impedance (default 1). Columns follow vsh_eval.m ([l m type], type 1 = TE, type 2 = TM), with
% A1 = r x grad Y / sqrt(l(l+1)), A2 = grad Y / sqrt(l(l+1)) on the unit sphere:
%   TE: E = z_l(kr) A1,           H = -(j/eta) Z_l(kr) A2
%   TM: E = Z_l(kr) A2,           H = -(j/eta) z_l(kr) A1
% with Z_l(x) = [x z_l(x)]' / x = z_{l-1}(x) - l z_l(x)/x. (From curl(f Phi) and curl(f Psi) of the
% vector spherical harmonics and H = (j/(k eta)) curl E; the TM sign is chosen so that E = +Z A2.)
% Far field of 'out': z_l -> j^(l+1) e^{-jkr}/(kr), Z_l -> j^l e^{-jkr}/(kr).
if nargin < 7, eta = 1; end
th = th(:); ph = ph(:); r = r(:) .* ones(size(th));
[A1t, A1p, modes] = vsh_eval(L, th, ph);         % type-1 columns hold A1 in (t, p); type 2 hold A2
x = k * r;
z = zeros(numel(x), L + 1); Z = z;                % z(:, l+1) = z_l(x)
for l = 0:L
    z(:, l+1) = sbes(l, x, kind);
end
for l = 1:L
    Z(:, l+1) = z(:, l) - l * z(:, l+1) ./ x;
end
nm = size(modes, 1);
Et = zeros(numel(x), nm); Ep = Et; Ht = Et; Hp = Et;
for q = 1:nm
    l = modes(q, 1); t = modes(q, 3);
    if t == 1                                     % TE: A1 in this column, A2 in the next
        a1t = A1t(:, q); a1p = A1p(:, q); a2t = A1t(:, q+1); a2p = A1p(:, q+1);
        Et(:, q) = z(:, l+1) .* a1t;  Ep(:, q) = z(:, l+1) .* a1p;
        Ht(:, q) = -1i/eta * Z(:, l+1) .* a2t;  Hp(:, q) = -1i/eta * Z(:, l+1) .* a2p;
    else                                          % TM: A2 in this column, A1 in the previous
        a2t = A1t(:, q); a2p = A1p(:, q); a1t = A1t(:, q-1); a1p = A1p(:, q-1);
        Et(:, q) = Z(:, l+1) .* a2t;  Ep(:, q) = Z(:, l+1) .* a2p;
        Ht(:, q) = -1i/eta * z(:, l+1) .* a1t;  Hp(:, q) = -1i/eta * z(:, l+1) .* a1p;
    end
end
end

function v = sbes(l, x, kind)
j = sqrt(pi ./ (2*x)) .* besselj(l + 0.5, x);
if strcmp(kind, 'reg'), v = j; return, end
y = sqrt(pi ./ (2*x)) .* bessely(l + 0.5, x);
v = j - 1i*y;                                     % h_l^(2): outgoing for e^{jwt}
end
