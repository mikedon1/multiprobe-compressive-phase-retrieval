function c = swe_fit(S, Et, Ep)
%SWE_FIT  Spherical-wave coefficients (vector, ordered as S.n/S.type/S.mi) of the tangential
%field Et (theta comp), Ep (phi comp), each nt x np on the grid of swe_setup.
np = numel(S.phi);
Ftm = (Et * conj(S.Ephi)) / np;                 % nt x (2N+1): phi Fourier components
Fpm = (Ep * conj(S.Ephi)) / np;
c = zeros(numel(S.n), 1);
for mi = 1:numel(S.m)
    c(S.mi == mi) = S.Pi{mi} * [Ftm(:, mi); Fpm(:, mi)];
end
end
