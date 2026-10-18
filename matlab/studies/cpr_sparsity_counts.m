function T = cpr_sparsity_counts()
%CPR_SPARSITY_COUNTS  Number of largest coefficients holding 99% of the energy of the example
%pattern (64 x 32 x 16) in several representations (Section II-D of the paper):
%   each    : a 2D DFT of each wavelength on its own (coefficients of all wavelengths pooled)
%   dft3    : a joint 3D DFT
%   dct     : 2D DFT x DCT along wavelength, no delay compensation
%   dctTau  : 2D DFT x DCT with the delay of Eq. (7) removed
%   rand3   : joint 3D DFT after a random phase on each wavelength (mean over 20 draws)
[X, f] = pattern_crlh(64, 16, 1);
[nw, nh, nd] = size(X);
k99 = @(c) find(cumsum(sort(abs(c(:)).^2, 'descend')) >= 0.99*sum(abs(c(:)).^2), 1);
F2 = zeros(size(X));
for k = 1:nd, F2(:, :, k) = fft2(X(:, :, k)) / sqrt(nw*nh); end
T.each = k99(F2);
T.dft3 = k99(fftn(X));
B = cpr_basis('dct', f, 0);                          % nd x nd unitary synthesis
proj = @(M) reshape(reshape(F2, [], nd) * conj(M), nw, nh, nd);   % s = M^H x along wavelength
T.dct = k99(proj(B));
x2 = reshape(X, [], nd);
tau = -angle(sum(sum(conj(x2(:, 1:end-1)) .* x2(:, 2:end)))) / (2*pi*(f(2) - f(1)));
T.tau_ns = tau*1e9;
T.dctTau = k99(proj(cpr_basis('dct', f, tau)));
rng(1); r = zeros(1, 20);
for i = 1:20, Xr = X .* reshape(exp(1i*2*pi*rand(1, nd)), 1, 1, nd); r(i) = k99(fftn(Xr)); end
T.rand3 = round(mean(r));
fprintf('99%% energy counts of %d coefficients: each %d, 3D DFT %d, DCT %d, DCT+delay %d (tau %.3f ns), random phases %d\n', ...
    numel(X), T.each, T.dft3, T.dct, T.dctTau, T.tau_ns, T.rand3);
end
