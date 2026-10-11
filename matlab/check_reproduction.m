function check_reproduction()
%CHECK_REPRODUCTION  Quick check that this installation reproduces a number of the paper.
%   >> setup; check_reproduction
% Runs the base phaseless IS configuration (paper preset, seed 1, about 30 s) and compares the
% NMSE with the value behind the point m/n = 0.3 of Fig. 3(c): -43.7 dB.
cpr_settings('paper');
r = cpr_run(cpr_cfg(), 1);
ref = -43.70; tol = 0.05;
ok = abs(r.score.global - ref) < tol;
fprintf('IS, m/n = 0.3, seed 1: NMSE %.2f dB (paper %.2f dB, tolerance %.2f dB): %s\n', ...
    r.score.global, ref, tol, ternary(ok, 'PASS', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end