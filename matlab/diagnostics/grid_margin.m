function grid_margin()
% Is the leg-free margin an AZIMUTH (rotation-step) condition or a general grid-oversampling condition?
% Rotation-phase margin (worst small per-rotation phase, fraction outside the span, dB) on unequal grids:
% nw theta samples over 0-360 (scan), nh rotations over 0-180, L = 14..20, truth = degree-16 fit, 4 wavelengths.
cd(fileparts(fileparts(fileparts(mfilename('fullpath'))))); setup; maxNumCompThreads(1); warning('off','all');
D = crlh_data(); nd = 4; Lfit = 16;
[T, Pp] = ndgrid(D.theta(:), D.phi(1:end-1)); keep = sin(T(:)) > 1e-6; w = sqrt(sin(T(keep)));
[Yt, Yp] = vsh_eval(Lfit, T(keep), Pp(keep)); A = [w .* Yt; w .* Yp];
B = zeros(2*nnz(keep), nd);
for t = 1:nd
    Et = D.E_theta{t}(:, 1:end-1); Ep = D.E_phi{t}(:, 1:end-1); B(:, t) = [w .* Et(keep); w .* Ep(keep)];
end
C = A \ B;
grids = [32 16; 32 24; 32 32; 48 16; 64 16; 48 24];
fprintf('%6s %6s %10s | margin (dB) at L = 14 16 18 20 22 | basis cols / 2 nw nh\n', 'nw', 'nh', 'azimuths');
for g = 1:size(grids, 1)
    nw = grids(g,1); nh = grids(g,2);
    [th, phs, sg] = scan(nw, nh);
    [Ys_t, Ys_p] = vsh_eval(Lfit, th, phs);
    Xs = [sg .* (Ys_t * C); sg .* (Ys_p * C)];                       % (2 n2) x nd
    n2 = nw*nh; r = (1:2*n2).'; col = floor((mod(r-1, n2))/nw) + 1;
    out = zeros(1, 5); cols = zeros(1, 5); Ls = [14 16 18 20 22];
    for i = 1:numel(Ls)
        Pb = basis(th, phs, sg, Ls(i)); cols(i) = size(Pb, 2);
        Mo = zeros(numel(Xs), nh); Mr = Mo;
        for l = 1:nh
            V = 1i*Xs .* (col == l); Mr(:, l) = V(:); W = V - Pb*(Pb'*V); Mo(:, l) = W(:);
        end
        Q = null(ones(1, nh)); Wo = [real(Mo); imag(Mo)]*Q; Wr = [real(Mr); imag(Mr)]*Q;
        [~, Rw] = qr(Wr, 0); s = sort(svd(Wo / Rw)); out(i) = 20*log10(s(1));
    end
    fprintf('%6d %6d %10d | %s | %s / %d\n', nw, nh, 2*nh, mat2str(round(out, 1)), mat2str(cols), 2*n2);
end
disp('GRID DONE')
end
function [th, phs, sg] = scan(nw, nh)
[ts, ps] = ndgrid((0:nw-1)*2*pi/nw, (0:nh-1)*pi/nh); ts = ts(:); ps = ps(:);
far = ts > pi + 1e-9; th = ts; th(far) = 2*pi - ts(far); p = ps; p(far) = ps(far) + pi;
phs = mod(p - pi/2, 2*pi); sg = ones(size(ts)); sg(far) = -1;
end
function Psi = basis(th, phs, sg, L)
[Yt, Yp] = vsh_eval(L, th, phs); Y = [sg .* Yt; sg .* Yp];
Psi = zeros(size(Y)); k = 0;
for j = 1:size(Y, 2)
    v = Y(:, j); n0 = norm(v); if n0 == 0, continue, end
    for pass = 1:2, v = v - Psi(:, 1:k) * (Psi(:, 1:k)' * v); end
    if norm(v) > 1e-6 * n0, k = k + 1; Psi(:, k) = v / norm(v); end
end
Psi = Psi(:, 1:k);
end
