function ok = cpr_test_operator()
%CPR_TEST_OPERATOR  Check SparseDftLinTrans against an explicitly formed dense operator.
%   ok = cpr_test_operator()
% Builds a small random sensing matrix Phi (sparse rows, complex unit-modulus weights such as
% calibrated cable phases and tone factors) and the dense H = Phi * kron(M, Psi2), where Psi2 is the
% unitary inverse 2D DFT and M the wavelength synthesis (DCT with a delay). It then compares
% mult, multTr, multSq (|H|.^2 v), and multSqTr (|H|'.^2 u) of the fast operator with the dense ones,
% for one and for two polarizations (M block diagonal). Errors above 1e-10 (relative) fail.
cpr_setup();
nw = 8; nh = 4; n2 = nw*nh;
rng(7);
cases = {1, 3; 2, 3};                                % {npol, nd per polarization}
ok = true;
for c = 1:size(cases, 1)
    npol = cases{c, 1}; nd1 = cases{c, 2}; nd = npol*nd1;
    f = 1.9e9 + 50e6*(0:nd1-1);
    M = cpr_basis('dct', f, 2.3e-9, npol);
    n = n2*nd; m = 40;
    Phi = sparse(m, n);
    for i = 1:m                                      % P = 4 probes, possibly at different wavelengths
        j = randperm(n, 4);
        Phi(i, j) = exp(1i*2*pi*rand(1, 4)); %#ok<SPRIX>
    end
    Psi2 = zeros(n2);
    for b = 1:n2, e = zeros(nw, nh); e(b) = 1; Psi2(:, b) = reshape(ifft2(e), [], 1) * sqrt(n2); end
    H = full(Phi) * kron(M, Psi2);
    A = SparseDftLinTrans(Phi, nw, nh, nd, M, 'exact');
    s = randn(n, 1) + 1i*randn(n, 1); z = randn(m, 1) + 1i*randn(m, 1);
    v = rand(n, 1); u = rand(m, 1);
    err = [rel(A.mult(s), H*s), rel(A.multTr(z), H'*z), ...
           rel(A.multSq(v), abs(H).^2*v), rel(A.multSqTr(u), (abs(H).^2)'*u)];
    pass = all(err < 1e-10);
    fprintf('npol %d: mult %.1e, multTr %.1e, multSq %.1e, multSqTr %.1e  %s\n', npol, err, tern(pass, 'PASS', 'FAIL'));
    ok = ok && pass;
end
if ~ok, error('cpr_test_operator:fail', 'SparseDftLinTrans does not match the dense operator.'); end
end

function e = rel(a, b), e = norm(a(:) - b(:)) / max(norm(b(:)), eps); end
function s = tern(c, a, b), if c, s = a; else, s = b; end, end
