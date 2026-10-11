function [Et, Ep] = swe_synth(S, c)
%SWE_SYNTH  Tangential field (theta and phi components, nt x np) on the swe_setup grid from
%spherical-wave coefficients c (vector, ordered as S.n/S.type/S.mi).
nt = numel(S.theta);
Ftm = zeros(nt, numel(S.m)); Fpm = Ftm;
for mi = 1:numel(S.m)
    v = S.A{mi} * c(S.mi == mi);
    Ftm(:, mi) = v(1:nt); Fpm(:, mi) = v(nt+1:end);
end
Et = Ftm * S.Ephi.'; Ep = Fpm * S.Ephi.';
end
