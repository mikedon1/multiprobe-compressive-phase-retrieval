function ok = cpr_test_vsh()
%CPR_TEST_VSH  Native correctness tests for the spherical-wave basis and the self-calibration step (2026-10-20).
%   ok = cpr_test_vsh()
% 1. Rotation law: rotating the field by psi about the polar axis, E'(th, ph) = E(th, ph - psi), multiplies the
%    physical mode coefficient c_lmt by exp(-i m psi) (vsh_eval.m modes).
% 2. Orthonormality of vsh_basis on the paper grid (32 x 16, both components, L = 12): Psi' * Psi = I.
% 3. Self-calibration: planted per-port cable lengths (zero mean) are recovered by cpr_selfcal from noise-free
%    amplitudes of a random multi-port operator, starting from zero, with no prior.
% Raises an error if any test fails; returns true otherwise.
rng(7);
% --- 1. rotation law ---------------------------------------------------------------------------------------
L = 6; psi = 0.37;
th = acos(2*rand(200, 1) - 1); ph = 2*pi*rand(200, 1);
[Yt, Yp, modes] = vsh_eval(L, th, ph);
[Yt2, Yp2] = vsh_eval(L, th, ph - psi);
c = randn(size(modes, 1), 1) + 1i*randn(size(modes, 1), 1);
cr = c .* exp(-1i*modes(:, 2)*psi);
e1 = norm([Yt2*c; Yp2*c] - [Yt*cr; Yp*cr]) / norm([Yt2*c; Yp2*c]);
fprintf('cpr_test_vsh: rotation law, relative error %.2e\n', e1);
% vsh_eval differentiates the Legendre functions numerically (step 1e-6), so allow 1e-8
assert(e1 < 1e-8, 'cpr_test_vsh: rotation law failed (%.2e)', e1);
% --- 2. orthonormal basis -------------------------------------------------------------------------------------
Psi = vsh_basis(32, 16, 12, 2);
e2 = norm(Psi'*Psi - eye(size(Psi, 2)), 'fro') / sqrt(size(Psi, 2));
fprintf('cpr_test_vsh: basis orthonormality, %d columns, error %.2e\n', size(Psi, 2), e2);
assert(size(Psi, 2) == 2*12*14, 'cpr_test_vsh: expected %d columns, got %d', 2*12*14, size(Psi, 2));
assert(e2 < 1e-10, 'cpr_test_vsh: basis not orthonormal (%.2e)', e2);
% --- 3. self-calibration recovers planted lengths ---------------------------------------------------------------
nP = 12; n = 60; m = 400; P = 4; nt = 3;
lam = [10.0; 10.9; 12.0];                              % cm
x = randn(n, 1) + 1i*randn(n, 1);
ii = []; jj = []; pp = []; tt = []; vv = [];
for r = 1:m
    ports = randperm(nP, P); tone = randi(nt);
    for q = 1:P
        j = randi(n);                                  % each port sees one pattern sample in this reading
        ii(end+1, 1) = r; jj(end+1, 1) = j; pp(end+1, 1) = ports(q); %#ok<AGROW>
        tt(end+1, 1) = tone; vv(end+1, 1) = exp(2i*pi*rand); %#ok<AGROW>
    end
end
ent = struct('i', ii, 'j', jj, 'port', pp, 'tone', tt, 'v', vv, 'nPorts', nP, 'm', m, 'n', n);
dTrue = 0.3*randn(nP, 1); dTrue = dTrue - mean(dTrue);    % cm, small enough for the Gauss-Newton basin
k = 1i*2*pi ./ lam(tt);
y = abs(accumarray(ii, vv .* exp(k .* dTrue(pp)) .* x(jj), [m 1]));
[~, d] = cpr_selfcal(ent, x, y, lam, zeros(nP, 1), 30, NaN);
e3 = norm((d - mean(d)) - dTrue) / norm(dTrue);
fprintf('cpr_test_vsh: self-calibration, relative length error %.2e\n', e3);
assert(e3 < 1e-6, 'cpr_test_vsh: self-calibration did not recover the lengths (%.2e)', e3);
ok = true;
end
