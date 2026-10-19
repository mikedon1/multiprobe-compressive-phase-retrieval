classdef GroupSparseScaEstim < SparseScaEstim
    %GROUPSPARSESCAESTIM  Bernoulli-Gaussian prior whose sparsity rate is learned per GROUP of
    %coefficients (structured sparsity), instead of one rate for all coefficients.
    %
    %   est = GroupSparseScaEstim(estim1, p1, x0, groups, mode, 'autoTune', true)
    %
    % groups: n x 1 positive integer labels. mode:
    %   'avg' : rate of coefficient j = mean posterior activity of its group (EM per group)
    %   'sep' : separable rate rho(j) = a(groups(j)) * b(cols(j)) / mean(activity), with a the group
    %           mean and b the mean over the coefficients sharing cols(j) (set obj.cols), e.g. a mode
    %           group times a frequency-basis index
    % The parent (SparseScaEstim, GAMPmatlab) computes the posterior activities py1 and its joint
    % rate; this class then replaces the joint rate by the structured one, which the next GAMP
    % iteration uses. Rates are clipped to [1e-6, 0.999].
    properties
        groups; mode = 'avg'; cols = [];
    end
    methods
        function obj = GroupSparseScaEstim(estim1, p1, x0, groups, mode, varargin)
            obj = obj@SparseScaEstim(estim1, p1, x0, varargin{:});
            obj.groups = groups(:); obj.mode = mode;
        end

        function [xhat, xvar, klDivNeg, py1] = estim(obj, rhat, rvar)
            tuning = obj.autoTune && ~obj.disableTune && obj.counter <= 0;
            [xhat, xvar, klDivNeg, py1] = estim@SparseScaEstim(obj, rhat, rvar);
            if ~tuning, return, end
            g = obj.groups; ng = max(g);
            a = accumarray(g, py1(:), [ng 1]) ./ accumarray(g, 1, [ng 1]);
            switch obj.mode
                case 'avg'
                    p = a(g);
                case 'sep'
                    c = obj.cols(:); nc = max(c);
                    b = accumarray(c, py1(:), [nc 1]) ./ accumarray(c, 1, [nc 1]);
                    p = a(g) .* b(c) / max(mean(py1(:)), 1e-12);
                otherwise
                    error('GroupSparseScaEstim: unknown mode %s', obj.mode);
            end
            obj.p1 = min(max(p, 1e-6), 0.999);
        end
    end
end
