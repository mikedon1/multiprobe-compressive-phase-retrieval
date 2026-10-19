classdef CalVshLinTrans < LinTrans
    %CALVSHLINTRANS  Section VI operator with embedded-probe calibration (GAMPmatlab LinTrans).
    %
    %   A = CalVshLinTrans(Phi, G, M)
    %
    % Unknowns: coefficients s (p x nf, frequency-basis along the columns) of the AUT's FREE-SPACE mode
    % coefficients a = s * M.' (p x nf: VSH mode q of sph_modes/vsh_eval order x frequency).
    % "Pattern" (what syn returns, scored and re-expressed by cpr_gauge): a(:), mode-major.
    % Probe responses: x_f = G{f} * a(:, f), with G{f}((l-1)*np + k, q) = R_kq(f) exp(-1i m_q psi_l), the
    % calibrated embedded response of port k to mode q at rotation l (fwc_probe_resp, fwc_predict), ordered
    % port fastest, then rotation, then frequency, the order of make_readings_fw. Readings: z = Phi * x.
    % mult/multTr: synthesis and its adjoint; ana: exact inverse of syn (M unitary); multSq/multSqTr exact
    % (|H|.^2 formed in single precision whenever M changes). Properties nd, n2 serve cpr_gauge.
    properties
        Phi; PhiH; G; GH; M; p; nf; nx; nd; n2; H2; H2T; L;
    end
    methods
        function obj = CalVshLinTrans(Phi, G, M)
            nf = numel(G); [nx, p] = size(G{1});
            obj = obj@LinTrans(size(Phi, 1), p*nf);
            if size(Phi, 2) ~= nx*nf, error('CalVshLinTrans: Phi has %d columns, expected %d', size(Phi, 2), nx*nf); end
            obj.Phi = Phi; obj.PhiH = Phi'; obj.G = G; obj.GH = cellfun(@(g) g', G, 'UniformOutput', false);
            obj.p = p; obj.nf = nf; obj.nx = nx; obj.nd = nf; obj.n2 = p;
            obj.L = round(sqrt(p/2 + 1) - 1);          % p = 2 L (L + 2) for the full dual basis (informative only)
            obj.setM(M);
        end

        function setM(obj, M)
            obj.M = M(1:obj.nf, 1:obj.nf);
            m = size(obj.Phi, 1); p_ = obj.p; nf_ = obj.nf; nx_ = obj.nx;
            PG = cell(1, nf_);                             % Phi_f * G_f (m x p), f = frequency
            for f = 1:nf_
                PG{f} = obj.Phi(:, (f-1)*nx_ + (1:nx_)) * obj.G{f};
            end
            obj.H2 = zeros(m, p_*nf_, 'single');
            for j = 1:nf_                                  % column block j = frequency-basis index
                Hj = 0;
                for f = 1:nf_
                    if obj.M(f, j) ~= 0, Hj = Hj + PG{f} * obj.M(f, j); end
                end
                obj.H2(:, (j-1)*p_ + (1:p_)) = single(abs(full(Hj)).^2);
            end
            obj.H2T = obj.H2.';
        end

        function a = syn(obj, s)
            a = reshape(reshape(s, obj.p, obj.nf) * obj.M.', [], 1);
        end

        function s = adj(obj, a)
            s = reshape(reshape(a, obj.p, obj.nf) * conj(obj.M), [], 1);
        end

        function s = ana(obj, a)
            s = obj.adj(a);                              % M unitary: the adjoint is the inverse
        end

        function x = probes(obj, a)
            % probe responses (nx*nf x 1) of mode coefficients a (p*nf x 1)
            A_ = reshape(a, obj.p, obj.nf); x = zeros(obj.nx, obj.nf);
            for f = 1:obj.nf, x(:, f) = obj.G{f} * A_(:, f); end
            x = x(:);
        end

        function z = mult(obj, s)
            z = obj.Phi * obj.probes(obj.syn(s));
        end

        function s = multTr(obj, z)
            X = reshape(obj.PhiH * z, obj.nx, obj.nf); A_ = zeros(obj.p, obj.nf);
            for f = 1:obj.nf, A_(:, f) = obj.GH{f} * X(:, f); end
            s = obj.adj(A_(:));
        end

        function q = multSq(obj, v)
            q = double(obj.H2 * single(v));
        end

        function s = multSqTr(obj, u)
            s = double(obj.H2T * single(u));
        end
    end
end
