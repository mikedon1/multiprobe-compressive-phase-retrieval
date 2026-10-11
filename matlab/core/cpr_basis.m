function M = cpr_basis(kind, f, tau, npol)
%CPR_BASIS  Unitary wavelength-axis synthesis matrix (nd x nd) with delay tau:
%   x_k = exp(-1i*2*pi*(f_k - f_1)*tau) * sum_f B(k,f) s_f
%   kind 'dct' : B = orthonormal DCT-II synthesis (C')
%   kind 'dft' : B = unitary inverse DFT
% A pure delay is a linear phase in frequency; removing it keeps a smooth
% pattern smooth along wavelength (the natural CRLH data step ~ -88 deg per
% 50 MHz, i.e. tau ~ 5 ns, which a DCT/DFT alone represents poorly).
% npol > 1 (dual polarization): block diagonal, the same M for each polarization's nd slices
% (one shared delay), i.e. kron(eye(npol), M).
if nargin < 3, tau = 0; end
if nargin < 4, npol = 1; end
nd = numel(f);
switch kind
    case 'dct'
        k = (0:nd-1).'; j = 0:nd-1;
        C = sqrt(2/nd) * cos(pi*(2*j+1).*k/(2*nd)); C(1,:) = C(1,:)/sqrt(2);
        B = C';
    case 'dft'
        B = fft(eye(nd))' / sqrt(nd);
    otherwise
        error('cpr_basis: unknown kind %s', kind);
end
M = diag(exp(-1i*2*pi*(f(:) - f(1))*tau)) * B;
if npol > 1, M = kron(eye(npol), M); end
end
