function ts_times()
%TS_TIMES  Two-scan solve times per frequency from the stored runs (random starts only, and with the continuation), totals including the
%unsuccessful random-start searches.

A = load(fullfile(cpr_paths().results, 'ift_floor_crlh_r15_30_snr50_randomstarts.mat')); B = load(fullfile(cpr_paths().results, 'ift_floor_crlh_r15_30_snr50.mat'));
Ra = A.R; Rb = B.R; fn = fieldnames(Ra); fprintf('fields: %s\n', strjoin(fn.', ' '));
[~, o] = sort([Ra.t]); Ra = Ra(o); [~, o] = sort([Rb.t]); Rb = Rb(o);
fprintf('tBuild %s | build field in B: %s\n', mat2str(A.tBuild), mat2str(B.tBuild));
fprintf(' t | random-only: starts  time(s)  res(dB) reached | final: init           starts  time(s)  res(dB) reached\n');
for k = 1:numel(Ra)
    a = Ra(k); b = Rb(k); ra = pick(a, {'res_dB','res','nres_dB','resdB'}); rb = pick(b, {'res_dB','res','nres_dB','resdB'});
    ia = 'random'; if isfield(b, 'init'), ib = b.init; else, ib = '?'; end
    fprintf('%2d | %6d %9.0f %8.1f %d | %-12s %6d %9.0f %8.1f\n', a.t, pick(a, {'starts','nstart'}), pick(a, {'time'}), ra, pick(a, {'reached'}), ib, pick(b, {'starts','nstart'}), pick(b, {'time'}), rb);
end
tA = [Ra.time]; tB = [Rb.time];
fprintf('TOTAL random-only (all 16 frequencies, including the unsuccessful 400-start searches): %.0f s = %.1f min = %.2f h\n', sum(tA), sum(tA)/60, sum(tA)/3600);
fprintf('TOTAL with continuation (final entries): %.0f s = %.1f min; the unsuccessful random searches that the continuation follows are in the random-only entries of those frequencies\n', sum(tB), sum(tB)/60);
end
function v = pick(s, names)
v = NaN; for i = 1:numel(names), if isfield(s, names{i}), v = s.(names{i}); if islogical(v), v = double(v); end, if numel(v) > 1, v = v(end); end, return, end, end
end