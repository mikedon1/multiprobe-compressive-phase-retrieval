function Psi = vsh_basis(nw, nh, L, npol)
%VSH_BASIS  Orthonormal angular basis of the measured components on the scan grid.
%   Psi = vsh_basis(nw, nh, L, npol)
% npol = 1 (default): the across component (sg * phi component, vsh_scan.m) of the far-field modes
%   up to degree L (vsh_eval.m), n2 x p. One tangential component spans about half the modes (for
%   each m, the phi components of both mode types lie in sin^(|m|-1)(theta) times polynomials in
%   cos(theta)), so the columns, ordered by degree, are orthonormalized by Gram-Schmidt (with
%   reorthogonalization) and the dependent ones skipped.
% npol = 2: both components, [along; across] stacked (2 n2 x p, along = sg * theta component, the
%   order of the dual-polarized arc ports in crlh_sphere.m / cpr_sphere_run.m). The
%   modes are then independent; the same orthonormalization gives p = 2L(L+2) columns.
% Psi has orthonormal columns on the scan grid; the leading columns are the low degrees.
if nargin < 4, npol = 1; end
persistent cache
key = sprintf('%d_%d_%d_%d', nw, nh, L, npol);
if isempty(cache), cache = containers.Map(); end
if isKey(cache, key), Psi = cache(key); return, end
[th, phs, sg] = vsh_scan(nw, nh);
[Yt, Yp] = vsh_eval(L, th, phs);                % columns ordered by degree l
Y = sg .* Yp;
if npol == 2, Y = [sg .* Yt; Y]; end
Psi = zeros(size(Y)); k = 0;
for j = 1:size(Y, 2)
    v = Y(:, j); n0 = norm(v);
    if n0 == 0, continue, end
    for pass = 1:2
        v = v - Psi(:, 1:k) * (Psi(:, 1:k)' * v);
    end
    if norm(v) > 1e-6 * n0
        k = k + 1; Psi(:, k) = v / norm(v);
    end
end
Psi = Psi(:, 1:k);
cache(key) = Psi;
end
