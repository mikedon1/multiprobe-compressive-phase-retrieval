function [c, Neff, spec] = swe_trunc(S, c, tol_dB)
%SWE_TRUNC  Keep only the mode degrees n whose energy is within tol_dB of the strongest degree
%(highest such n = Neff); zero the rest. Needed before propagating INWARD (far -> near, or
%r2 -> r1): the radial ratio grows like (2n-1)!!/(kr)^(n+1) for n > kr, so numerical noise in
%degrees the antenna does not radiate would otherwise be amplified by up to ~1e12. Neff is
%also the minimum-sphere limit (Neff ~ k r0 + a few), i.e. the spherical aperture constraint.
% spec: per-degree energy (dB re the strongest degree), 1..S.N.
spec = accumarray(S.n, abs(c).^2, [S.N 1]);
spec = 10*log10(spec / max(spec));
Neff = find(spec > tol_dB, 1, 'last');
c(S.n > Neff) = 0;
end
