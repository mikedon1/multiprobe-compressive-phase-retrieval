function cpr_robustness_report()
%CPR_ROBUSTNESS_REPORT  Print the numbers of the paper's Gaussian-weight comparison (Section IV-A) and
%cable-calibration-error sensitivity (Section IV-C) from the studies pr_gauss and pr_calerr.
%   run_study('pr_gauss'); run_study('pr_calerr'); cpr_robustness_report
% NMSE is the whole-pattern value with one global phase removed (dB).
T = cpr_load('pr_gauss');
fprintf('\nOn/off switches vs complex Gaussian weights (NMSE in dB, one trial each)\n');
fprintf('%-6s %8s %10s %10s\n', 'arch', 'm/n', 'switches', 'Gaussian');
for a = {'IS', 'IC'}
    s = T(strcmp(T.arch, a{1}), :);
    for rt = unique(s.ratio).'
        sw = s.global(s.ratio == rt & strcmp(s.weights, 'switch'));
        ga = s.global(s.ratio == rt & strcmp(s.weights, 'gauss'));
        if ~isempty(sw) && ~isempty(ga), fprintf('%-6s %8.2f %10.2f %10.2f\n', a{1}, rt, sw(1), ga(1)); end
    end
end
C = cpr_load('pr_calerr');
fprintf('\nSensitivity to rms cable phase errors (NMSE in dB, one trial each)\n');
fprintf('%-6s %10s %10s %14s\n', 'arch', 'error (deg)', 'NMSE', '20log10(sigma)');
for a = {'IS', 'IC'}
    s = sortrows(C(strcmp(C.arch, a{1}), :), 'calErrDeg');
    for i = 1:height(s)
        lv = NaN; if s.calErrDeg(i) > 0, lv = 20*log10(s.calErrDeg(i)*pi/180); end
        fprintf('%-6s %10g %10.2f %14.1f\n', a{1}, s.calErrDeg(i), s.global(i), lv);
    end
end
end