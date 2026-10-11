function g = swe_radial(S, k, r)
%SWE_RADIAL  Radial factor of every coefficient of S at radius r (exp(+j w t), outgoing waves):
%   TE: h_n^(2)(kr)          TM: (1/kr) d/d(kr) [kr h_n^(2)(kr)] = h_{n-1}^(2)(kr) - n h_n^(2)(kr)/(kr)
% r = Inf gives the far-field form with the exp(-j k r)/(k r) factor removed and r = 1 m
% put back, i.e. the coefficients of a far field E(theta,phi) given as the asymptotic field
% evaluated at r = 1 m (the openEMS nf2ff convention, exp(-j k r)/r included):
%   TE: 1j^(n+1) exp(-j k)/k        TM: 1j^n exp(-j k)/k
% Coefficients at radius r2 from radius r1:  c2 = c1 .* swe_radial(S,k,r2) ./ swe_radial(S,k,r1).
n = S.n; te = S.type == 1;
if isinf(r)
    g = (1j.^(n + te)) * exp(-1j*k) / k;
    return
end
x = k*r;
hn  = sph_h2(n, x); hn1 = sph_h2(n - 1, x);
g = hn;
g(~te) = hn1(~te) - n(~te) .* hn(~te) / x;
end

function h = sph_h2(n, x)
h = sqrt(pi/(2*x)) * besselh(n + 0.5, 2, x);
end
