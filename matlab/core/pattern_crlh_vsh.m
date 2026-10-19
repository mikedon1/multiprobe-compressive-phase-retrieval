function [X, f] = pattern_crlh_vsh(nw, nd, f0idx, Lfit, npol, src)
%PATTERN_CRLH_VSH  The CRLH pattern of pattern_crlh.m, evaluated exactly at the scan points.
%   [X, f] = pattern_crlh_vsh(nw, nd, f0idx, Lfit, npol)
% The openEMS far field (both components, 2-deg grid) is fitted with vector spherical harmonics
% up to degree Lfit (default 16; the fit residual is below -140 dB) and evaluated at the scan
% directions (vsh_scan.m), instead of linearly interpolating the 2-deg grid.
%   npol = 1 (default): the across component, nw x nh x nd
%   npol = 2          : [along, across] stacked along dimension 3 (nw x nh x 2nd: along f_1..f_nd,
%                       then across f_1..f_nd), the order of vsh_basis(..., 2); along = theta
%                       component, reversed on the far half
% Same normalization as pattern_crlh (max |X| = 1 over everything returned).
if nargin < 3, f0idx = 1; end
if nargin < 4, Lfit = 16; end
if nargin < 5, npol = 1; end
nh = nw/2;
% src (2026-10-19): 'crlh' (default) | 'array' (array_data, the electrically larger second example) | a struct
% in the crlh_data format
if nargin < 6 || isempty(src), src = 'crlh'; end
if isstruct(src), D = src;
elseif startsWith(src, 'array')                          % 'array' = 5 x 5, 'array4' = 4 x 4, ...
    na = str2double(src(6:end)); if isnan(na), na = 5; end
    D = array_data(na, na, 6, 20);
else, D = crlh_data(); end
[T, Pp] = ndgrid(D.theta(:), D.phi(1:end-1));            % source grid without the 360-deg column
keep = sin(T(:)) > 1e-6;
w = sqrt(sin(T(keep)));
[Yt, Yp] = vsh_eval(Lfit, T(keep), Pp(keep));
A = [w .* Yt; w .* Yp];
idx = f0idx + (0:nd-1);
B = zeros(2*nnz(keep), nd);
for t = 1:nd
    Et = D.E_theta{idx(t)}(:, 1:end-1); Ep = D.E_phi{idx(t)}(:, 1:end-1);
    B(:, t) = [w .* Et(keep); w .* Ep(keep)];
end
C = A \ B;
res = norm(A*C - B, 'fro') / norm(B, 'fro');
if res > 1e-5, warning('pattern_crlh_vsh: fit residual %.1f dB', 20*log10(res)); end
[th, phs, sg] = vsh_scan(nw, nh);
[Ys_t, Ys_p] = vsh_eval(Lfit, th, phs);
X = reshape(sg .* (Ys_p * C), nw, nh, nd);
if npol == 2
    X = cat(3, reshape(sg .* (Ys_t * C), nw, nh, nd), X);
end
X = X / max(abs(X(:)));
f = D.freq(idx); f = f(:).';
end
