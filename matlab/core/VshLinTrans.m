classdef VshLinTrans < LinTrans
    %VSHLINTRANS  H = Phi * Psi with a vector-spherical-harmonic angular basis (GAMPmatlab LinTrans;
    %same interface as SparseDftLinTrans).
    %
    %   Phi : sparse m x (n2*nd) pattern-domain sensing matrix (make_readings.m); nd = npol*nf
    %         slices ordered [pol 1: f_1..f_nf, pol 2: f_1..f_nf]
    %   Psi : (npol*n2) x p angular synthesis on the scan grid (vsh_basis.m), orthonormal columns;
    %         both polarizations share the mode coefficients
    %   M   : nf x nf unitary wavelength-axis synthesis with delay (cpr_basis.m; a block-diagonal
    %         kron(eye(npol), M) is accepted and its first block used)
    %
    % mult/multTr are the synthesis and its ADJOINT; ana is the least-squares INVERSE of syn
    % (cpr_gauge re-expresses estimates with it). multSq/multSqTr are exact: |H|.^2 is formed
    % (single precision) whenever M changes.
    properties
        Phi; PhiH; Psi; PsiPinv; nw; nh; nd; nf; npol; n2; p; M; H2; H2T; L;
    end
    methods
        function obj = VshLinTrans(Phi, nw, nh, nd, M, L, npol, Psi)
            % Psi (optional): a basis for another grid, e.g. the full-wave arc (vsh_basis_grid.m)
            if nargin < 7, npol = 1; end
            if nargin < 8 || isempty(Psi), Psi = vsh_basis(nw, nh, L, npol); end
            nf = nd / npol;
            obj = obj@LinTrans(size(Phi, 1), size(Psi, 2)*nf);
            obj.Phi = Phi; obj.PhiH = Phi'; obj.Psi = Psi; obj.PsiPinv = pinv(Psi);
            obj.nw = nw; obj.nh = nh; obj.nd = nd; obj.nf = nf; obj.npol = npol;
            obj.n2 = nw*nh; obj.p = size(Psi, 2); obj.L = L;
            obj.setM(M);
        end

        function setM(obj, M)
            nf_ = obj.nf; n2_ = obj.n2; np = obj.p;
            M = M(1:nf_, 1:nf_);
            obj.M = M;
            obj.H2 = zeros(size(obj.Phi, 1), np*nf_, 'single');
            I2 = speye(n2_);
            for f = 1:nf_
                Hf = 0;
                for q = 1:obj.npol
                    Phiq = obj.Phi(:, (q-1)*nf_*n2_ + (1:nf_*n2_));
                    Hf = Hf + (Phiq * kron(M(:, f), I2)) * obj.Psi((q-1)*n2_ + (1:n2_), :);
                end
                obj.H2(:, (f-1)*np + (1:np)) = single(abs(full(Hf)).^2);
            end
            obj.H2T = obj.H2.';
        end

        function setPhi(obj, Phi)
            % new pattern-domain matrix (self-calibration); |H|.^2 is recomputed
            obj.Phi = Phi; obj.PhiH = Phi'; obj.setM(obj.M);
        end

        function x = syn(obj, s)
            Y = obj.Psi * reshape(s, obj.p, obj.nf) * obj.M.';          % (npol n2) x nf
            x = reshape(permute(reshape(Y, obj.n2, obj.npol, obj.nf), [1 3 2]), [], 1);
        end

        function s = adj(obj, x)
            s = reshape(obj.Psi' * obj.stack(x) * conj(obj.M), [], 1);
        end

        function s = ana(obj, x)
            s = reshape(obj.PsiPinv * obj.stack(x) * conj(obj.M), [], 1);
        end

        function Y = stack(obj, x)
            % pattern vector (slices [pol 1 f_1..f_nf, pol 2 ...]) -> (npol n2) x nf
            Y = reshape(permute(reshape(x, obj.n2, obj.nf, obj.npol), [1 3 2]), [], obj.nf);
        end

        function z = mult(obj, s)
            z = obj.Phi * obj.syn(s);
        end

        function s = multTr(obj, z)
            s = obj.adj(obj.PhiH * z);
        end

        function q = multSq(obj, v)
            q = double(obj.H2 * single(v));
        end

        function s = multSqTr(obj, u)
            s = double(obj.H2T * single(u));
        end
    end
end
