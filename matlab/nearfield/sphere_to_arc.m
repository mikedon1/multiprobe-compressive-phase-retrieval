function [Xt, Xp] = sphere_to_arc(Et, Ep, map, arcSign)
%SPHERE_TO_ARC  Tangential field on the theta x phi tensor grid -> the two components measured
%by dual-polarized probes on the paper's great-circle grid (arc tangent and arc normal), with
%map = D.map and arcSign = D.arcSign from crlh_sphere.m / array_sphere.m. Used for both the
%true patterns and the two-scan estimates, so both methods are scored on the same quantities.
Xt = arcSign .* Et(map);
Xp = arcSign .* Ep(map);
end
