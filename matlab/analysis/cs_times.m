function cs_times()
%CS_TIMES  Times stored with the one-sphere compressive runs of the two-scan comparison.

d = dir(fullfile(cpr_paths().results, 'cpr_sphere_r15_s1_snr50_pol2*.mat'));
for i = 1:numel(d)
    S = load(fullfile(d(i).folder, d(i).name)); fn = fieldnames(S); fprintf('%s: %s\n', d(i).name, strjoin(fn.', ' '));
    for j = 1:numel(fn)
        v = S.(fn{j});
        if isstruct(v)
            f2 = fieldnames(v); fprintf('   %s fields: %s\n', fn{j}, strjoin(f2.', ' '));
            for q = {'time', 'timeSetup', 'seconds'}, if isfield(v, q{1}), fprintf('   %s.%s = %s\n', fn{j}, q{1}, mat2str(double(v(1).(q{1})(:).'))); end, end
        end
    end
end
end