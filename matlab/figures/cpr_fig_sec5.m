function cpr_fig_sec5()
%CPR_FIG_SEC5  Section V figures in the paper's format (cpr_fig: 12 x 6 cm, Times, LaTeX axis labels).
%  ift_freq.pdf  : NMSE per frequency at 50 dB SNR, two-sphere IFT vs. dual-pol CPR (IC). Both are
%                  scored on both tangential components on the first sphere, one phase per
%                  wavelength (IFT: R(t).fieldErr; CPR: r.scoreVec.each).
%  ift_accel.pdf : acceleration t_IFT/t_CPR vs. gamma = t_m/t_s from eqs. (mtime_ift_norm) and
%                  (mtime_cpr), at n_lambda = 16, K = 4, rho = 0.6 (the Section V values), for
%                  several grid sizes (n_theta = n_phi/2).
outDir = cpr_paths().results;
pDir = fullfile(outDir, 'figs', 'paper');
F = cpr_fig();
% --- NMSE per frequency ---
S = load(fullfile(outDir, 'ift_floor_crlh_r15_30_snr50.mat')); R = S.R;
[~, o] = sort([R.t]); R = R(o);
fG = [R.f]/1e9; eI = [R.fieldErr];
% CPR estimates are mapped to the two-scan tensor grid (cpr_to_tensor) and scored exactly like
% R.fieldErr: both tangential components, one phase per frequency (all on one grid).
D1 = crlh_sphere(0.15, cpr_settings().nw, cpr_settings().nd, 1);
runs = {'cpr_sphere_r15_s1_snr50_pol2.mat', 'IC, {\itm}/{\itn} = 0.6'; ...
        'cpr_sphere_r15_s1_snr50_pol2_IS_m60.mat', 'IS, {\itm}/{\itn} = 0.6'; ...
        'cpr_sphere_r15_s1_snr50_pol2_IS_m30.mat', 'IS, {\itm}/{\itn} = 0.3'};
eC = zeros(size(runs, 1), numel(fG));
for i = 1:size(runs, 1)
    C = load(fullfile(outDir, runs{i, 1})); r = C.r;
    [Et, Ep] = cpr_to_tensor(double(r.Xh), D1);
    for t = 1:numel(fG)
        a = [reshape(Et(:,:,t), [], 1); reshape(Ep(:,:,t), [], 1)];
        b = [reshape(D1.Et(:,:,t), [], 1); reshape(D1.Ep(:,:,t), [], 1)];
        a = a * exp(1i*angle(a'*b)); eC(i, t) = 20*log10(norm(a - b)/norm(b));
    end
    a = []; b = [];
    for t = 1:numel(fG)                                   % whole-pattern NMSE, same alignment
        at = [reshape(Et(:,:,t), [], 1); reshape(Ep(:,:,t), [], 1)];
        bt = [reshape(D1.Et(:,:,t), [], 1); reshape(D1.Ep(:,:,t), [], 1)];
        a = [a; at*exp(1i*angle(at'*bt))]; b = [b; bt]; %#ok<AGROW>
    end
    fprintf('%-14s m = %d: per frequency %s, whole pattern %.2f dB\n', runs{i, 2}, r.m, ...
        mat2str(round(eC(i, :), 1)), 20*log10(norm(a - b)/norm(b)));
end
f = F.new(); ax = axes(f); hold(ax, 'on');
plot(ax, fG, eI, F.lss{1}, 'Color', F.cols(1,:), 'LineWidth', F.lw);
for i = 1:size(runs, 1)
    plot(ax, fG, eC(i, :), F.lss{i+1}, 'Color', F.cols(i+1,:), 'LineWidth', F.lw);
end
xlabel(ax, 'Frequency (GHz)'); ylabel(ax, 'NMSE (dB)');
legend(ax, [{'Two-scan'}, runs(:, 2).'], 'Location', 'northwest');
xlim(ax, [fG(1) fG(end)]); F.style(ax); F.save(f, fullfile(pDir, 'ift_freq'));
fprintf('two-scan %s\n', mat2str(round(eI, 1)));
% --- acceleration vs gamma ---
nl = 16; K = 4; rho = 0.6; g = logspace(0, 4, 200); nphis = [64 48 32 16];
f = F.new(); ax = axes(f); hold(ax, 'on'); nm = {};
for i = 1:numel(nphis)
    np = nphis(i); nt = 2*np;    % arc probes x rotations as in the IC layout of Fig. 13 (64 x 32)
    tI = g*(2*np - 1) + 4*nt*np*nl;
    tC = g*(np - 1) + 2*rho*nt*np*nl/K;
    plot(ax, g, tI./tC, F.lss{i}, 'Color', F.cols(i,:), 'LineWidth', F.lw);
    nm{end+1} = sprintf('{\\itn}_\\theta \\times {\\itn}_\\phi = %d \\times %d', nt, np); %#ok<AGROW>
end
set(ax, 'XScale', 'log'); xlim(ax, g([1 end]));
xlabel(ax, '$\gamma$', 'Interpreter', 'latex'); ylabel(ax, 'CPR acceleration');
legend(ax, nm, 'Location', 'southwest');   % plain TeX: LaTeX legends are sized too narrow (clipped)
F.style(ax); F.save(f, fullfile(pDir, 'ift_accel'));
end
