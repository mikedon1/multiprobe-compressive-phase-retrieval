function ok = check_reproduction()
%CHECK_REPRODUCTION  Quick check of this installation (about 1 minute).
%   >> setup; check_reproduction
% 1. Invariant test: the fast sensing operator matches an explicitly formed dense operator
%    (cpr_test_operator).
% 2. Regression test: the base phaseless IS configuration (paper preset, seed 1) reproduces the
%    NMSE behind the point m/n = 0.3 of the paper's sampling-ratio figure.
% Raises an error if either check fails.
cpr_test_operator();
cpr_settings('paper');
r = cpr_run(cpr_cfg(), 1);
ref = -49.96; tol = 0.05;
ok = abs(r.score.global - ref) < tol;
fprintf('IS, m/n = 0.3, seed 1: NMSE %.2f dB (reference %.2f dB, tolerance %.2f dB): %s\n', ...
    r.score.global, ref, tol, ternary(ok, 'PASS', 'FAIL'));
if ~ok
    error('check_reproduction:fail', ['The regression value differs from the reference. Results ' ...
        'can differ slightly across MATLAB versions and numerical libraries.']);
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
