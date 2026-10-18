function [Et, Ep] = cpr_to_tensor(Xh, D)
%CPR_TO_TENSOR  Map a dual-polarized CPR estimate on the paper's 64 x 32 great-circle grid (arc
%components, Xh = nw x nh x 2*nd: [theta-port slices, phi-port slices]) to the tangential field
%E_theta, E_phi on the (nw/2+1) x nw theta x phi tensor grid of crlh_sphere.m, the grid the
%two-scan technique is scored on (Section VII, Fig. 10). Inverse of sphere_to_arc.m.
%Away from the poles the two grids share every point (D.map is one-to-one). Each pole is hit at
%nh of its nw phi samples; the tensor sample at phi + 180 deg is the same field vector with both
%unit vectors reversed, so it is filled as minus the sample at phi.
[nw, nh, ns] = size(Xh); nd = ns/2;
nt = nw/2 + 1;
Et = complex(nan(nt, nw, nd)); Ep = Et;
for t = 1:nd
    a = complex(nan(nt, nw)); b = a;
    a(D.map) = D.arcSign .* Xh(:,:,t);
    b(D.map) = D.arcSign .* Xh(:,:,nd + t);
    for r = [1 nt]                                   % pole rows
        for j = 1:nw
            if isnan(a(r, j))
                jo = mod(j - 1 + nw/2, nw) + 1;      % phi + 180 deg
                a(r, j) = -a(r, jo); b(r, j) = -b(r, jo);
            end
        end
    end
    Et(:,:,t) = a; Ep(:,:,t) = b;
end
if any(isnan(Et(:))) || any(isnan(Ep(:))), error('cpr_to_tensor: unfilled tensor points'); end
end
