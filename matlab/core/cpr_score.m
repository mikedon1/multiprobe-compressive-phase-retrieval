function s = cpr_score(X, xh, nd)
%CPR_SCORE  NMSE of a pattern estimate (dB), three ways:
%   raw       : no phase alignment (the right metric for complex measurements,
%               which carry an absolute phase reference)
%   global    : one global phase removed
%   perlambda : one phase per wavelength removed (paper eq. align2; the right
%               metric when per-wavelength phases are unobservable, K = 1)
% Also per-wavelength NMSE (perlambda alignment). X: n2 x nd or nw x nh x nd.

X = reshape(X, [], nd);
if isempty(xh)
    s = struct('raw', NaN, 'global', NaN, 'perlambda', NaN, 'each', NaN(1, nd)); return
end
Xh = reshape(xh, [], nd);
nr = @(E) 20*log10(norm(E(:))/norm(X(:)));
s.raw = nr(X - Xh);
s.global = nr(X - Xh*exp(1i*angle(Xh(:)'*X(:))));
s.each = zeros(1, nd);
for k = 1:nd
    Xh(:,k) = Xh(:,k)*exp(1i*angle(Xh(:,k)'*X(:,k)));
    s.each(k) = 20*log10(norm(X(:,k) - Xh(:,k))/norm(X(:,k)));
end
s.perlambda = nr(X - Xh);
end
