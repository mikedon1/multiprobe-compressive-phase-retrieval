function [Psi, nraw] = vsh_basis_grid(th, ph, L, comp, keepMode)
%VSH_BASIS_GRID  Orthonormal VSH basis of one field component on an arbitrary set of directions.
%   [Psi, nraw] = vsh_basis_grid(th, ph, L, comp, keepMode)
% th, ph (rad): sample directions (any order; rows of Psi follow it). comp: 'theta' or 'phi'.
% keepMode (optional): @(modes) logical selection of modes ([l m type] rows, vsh_eval.m), e.g. one
% inversion parity. The component of the selected modes up to degree L is orthonormalized on the
% samples by Gram-Schmidt in degree order (with reorthogonalization); columns whose new part is
% below 1e-6 of their norm are skipped, as for vsh_basis.m. nraw: number of selected modes.
[Yt, Yp, modes] = vsh_eval(L, th, ph);
if strcmp(comp, 'theta'), Y = Yt; else, Y = Yp; end
if nargin >= 5 && ~isempty(keepMode), Y = Y(:, keepMode(modes)); end
nraw = size(Y, 2);
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
end
