function e = cpr_gauge(e, A, mode, f, kind, rowTone)
%CPR_GAUGE  Re-choose the representation of the current GAMP estimate between
%EM rounds. Neither mode changes the data fit; both only make the
%coefficients sparser.
%
%   mode 'perlambda' (phaseless, K = 1): each wavelength's global phase is
%       unobservable, so rotate wavelength k by exp(-1i*a_k) to align it with
%       wavelength k-1 (smoothest sequence). The pattern estimate changes by
%       the unobservable per-wavelength phases only; the measurement-domain
%       dual (shatNext) is rotated with its row's tone (rowTone).
%   mode 'delay' (complex, or phaseless K > 1): the relative phases are
%       measured, so the pattern is NOT changed. Instead the basis delay tau
%       is re-estimated from the estimate's mean phase step between adjacent
%       wavelengths, A's synthesis matrix is updated (A is a handle), and the
%       coefficients are re-expressed in the new basis.
nd = A.nd; n2 = A.n2;
if nd < 2, return; end
X = reshape(A.syn(e.xhat), n2, nd);
switch mode
    case 'perlambda'
        if nd ~= numel(f), error('cpr_gauge: perlambda (K = 1) not implemented for dual polarization'); end
        a = zeros(1, nd);
        for k = 2:nd
            a(k) = a(k-1) + angle(X(:,k-1)' * X(:,k));
        end
        d = kron(exp(-1i*a(:)), ones(n2, 1));
        rot = @(s) A.ana(d .* A.syn(s));
        e.xhat = rot(e.xhat); e.xhatNext = rot(e.xhatNext); e.xhatDamp = rot(e.xhatDamp);
        a = a(:);
        e.shatNext = exp(-1i*a(rowTone(:))) .* e.shatNext;
    case 'delay'
        % dual polarization: A.nd = npol*numel(f) slices; steps only within each polarization
        npol = round(nd / numel(f)); nf = numel(f); acc = 0;
        for q = 1:npol
            Xq = X(:, (q-1)*nf + (1:nf));
            acc = acc + sum(sum(conj(Xq(:,1:end-1)) .* Xq(:,2:end)));
        end
        step = angle(acc);                                            % mean phase step per tone
        tau = -step / (2*pi*(f(2) - f(1)));
        xs = {A.syn(e.xhat), A.syn(e.xhatNext), A.syn(e.xhatDamp)};
        A.setM(cpr_basis(kind, f, tau, npol));
        e.xhat = A.ana(xs{1}); e.xhatNext = A.ana(xs{2}); e.xhatDamp = A.ana(xs{3});
        e.tau = tau;
    otherwise
        error('cpr_gauge: unknown mode %s', mode);
end
end
