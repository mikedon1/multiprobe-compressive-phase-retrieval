classdef SparseDftLinTrans < LinTrans
    %SPARSEDFTLINTRANS  H = Phi * Psi without forming H (GAMPmatlab LinTrans).
    %
    %   Phi : sparse m x n pattern-domain sensing matrix (make_readings.m),
    %         columns a + (k-1)*n2 (angle a, wavelength k)
    %   Psi : unitary synthesis  x(:,k) = sum_f M(k,f) * Psi2 * s(:,f),
    %         Psi2 = unitary 2D inverse DFT over the nw x nh angle grid,
    %         M    = nd x nd unitary wavelength-axis synthesis (DCT' or
    %                inverse DFT, times the delay diag(exp(-1i*2*pi*(f-f1)*tau)))
    %
    % mult/multTr cost FFTs plus one sparse product. multSq/multSqTr are EXACT
    % (sqMode 'exact') using
    %   Psi2(a,b)*conj(Psi2(a',b)) = exp(+1i*2*pi*<a-a',b>)/n2,
    % so |H|.^2 * v needs, per row, only the (nnz per row)^2 pairs of its
    % entries plus one 2D FFT per wavelength pair (verified against the dense
    % H by cpr_test_operator.m). sqMode 'uniform' replaces |H_ij|^2 by its row
    % average ||Phi_i||^2/n (GAMP's uniform-variance approximation; faster).
    properties
        Phi; PhiH; nw; nh; nd; n2; M; Q; sqMode;
        Pm; PmT; rowE;   % Pm: m x (n2*nU) sparse pair matrix for exact |H|^2 products
        tUsed = [];      % wavelength pairs (k1,k2) that occur in some row (rows of Q used)
    end
    methods
        function obj = SparseDftLinTrans(Phi, nw, nh, nd, M, sqMode)
            [m, n] = size(Phi);
            obj = obj@LinTrans(m, n);
            if nargin < 6, sqMode = 'exact'; end
            obj.Phi = Phi; obj.PhiH = Phi'; obj.nw = nw; obj.nh = nh; obj.nd = nd;
            obj.n2 = nw*nh; obj.sqMode = sqMode;
            obj.rowE = full(sum(abs(Phi).^2, 2));
            if strcmp(sqMode, 'exact'), obj.buildPairs(); end
            obj.setM(M);
        end

        function setM(obj, M)
            % wavelength-axis synthesis matrix (nd x nd unitary)
            obj.M = M;
            nd_ = obj.nd;
            [k1, k2] = ndgrid(1:nd_, 1:nd_);                  % pair index t = k1 + (k2-1)*nd
            obj.Q = M(k1(:), :) .* conj(M(k2(:), :));           % nd^2 x nd
        end

        function buildPairs(obj)
            [ci, ri, vi] = find(obj.Phi.');                     % grouped by original row
            m = size(obj.Phi, 1);
            cnt = accumarray(ri, 1, [m 1]);
            start = cumsum([1; cnt(1:end-1)]);
            ce = cnt(ri);
            J1 = repelem((1:numel(ri)).', ce);
            within = (1:numel(J1)).' - repelem(cumsum([0; ce(1:end-1)]), ce);
            J2 = repelem(start(ri), ce) + within - 1;
            % For real v (variances), pair (J2,J1) contributes the complex
            % conjugate of pair (J1,J2), so keep J1 <= J2, double J1 < J2, and
            % take the real part: still exact, half the work.
            keep = J1 <= J2;
            wgt = 1 + (J1(keep) < J2(keep));
            J1 = J1(keep); J2 = J2(keep);
            n2_ = obj.n2; nw_ = obj.nw; nh_ = obj.nh; nd_ = obj.nd;
            a = mod(ci - 1, n2_); k = floor((ci - 1)/n2_);        % 0-based angle, tone
            it = mod(a, nw_); jp = floor(a/nw_);
            di = mod(it(J1) - it(J2), nw_); dj = mod(jp(J1) - jp(J2), nh_);
            t12 = k(J1) + k(J2)*nd_;                              % 0-based pair index
            [obj.tUsed, ~, tc] = unique(t12);                   % only wavelength pairs that occur
            obj.tUsed = obj.tUsed + 1;                            % 1-based rows of Q
            pIdx_ = di + dj*nw_ + (tc - 1)*n2_ + 1;
            % pair list as a sparse matrix (duplicates summed): one sparse
            % mat-vec per |H|^2 product instead of gather + accumarray
            obj.Pm = sparse(ri(J1), pIdx_, wgt .* vi(J1) .* conj(vi(J2)), m, n2_*numel(obj.tUsed));
            obj.PmT = obj.Pm.';                                   % non-conjugate transpose
        end

        function x = syn(obj, s)
            S = reshape(s, obj.nw, obj.nh, obj.nd);
            X = reshape(ifft2(S) * sqrt(obj.n2), obj.n2, obj.nd) * obj.M.';
            x = X(:);
        end

        function s = ana(obj, x)
            X = reshape(x, obj.n2, obj.nd) * conj(obj.M);
            s = reshape(fft2(reshape(X, obj.nw, obj.nh, obj.nd)) / sqrt(obj.n2), [], 1);
        end

        function z = mult(obj, s)
            z = obj.Phi * obj.syn(s);
        end

        function s = multTr(obj, z)
            s = obj.ana(obj.PhiH * z);
        end

        function p = multSq(obj, v)
            if strcmp(obj.sqMode, 'uniform')
                p = obj.rowE * (sum(v)/numel(v)); return
            end
            W = reshape(v, obj.n2, obj.nd) * obj.Q(obj.tUsed, :).'; % n2 x nU
            G = ifft2(reshape(W, obj.nw, obj.nh, []));
            p = real(obj.Pm * G(:));
        end

        function s = multSqTr(obj, u)
            n = obj.n2*obj.nd;
            if strcmp(obj.sqMode, 'uniform')
                s = ones(n, 1) * (sum(obj.rowE .* u)/n); return
            end
            C = obj.PmT * u;
            E = reshape(ifft2(reshape(full(C), obj.nw, obj.nh, [])), obj.n2, []);
            s = real(reshape(E * obj.Q(obj.tUsed, :), [], 1));
        end
    end
end
