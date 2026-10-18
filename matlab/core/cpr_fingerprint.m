function fp = cpr_fingerprint()
%CPR_FINGERPRINT  Settings and data identity stored with every checkpoint.
% A checkpoint is only resumed if its fingerprint equals the current one: the settings that change
% results (not paths, threads, or figure options), the pattern generator version, and the size of the
% data file. Raise dataVersion whenever the pattern generator or the data change.
S = cpr_settings();
S = rmfield(S, intersect(fieldnames(S), {'resultsDir', 'threads', 'figures', 'gampPath'}));
fp.settings = S;
fp.dataVersion = 2;              % 2: pattern_crlh shifts phi on the unique 180 columns (seam fix)
d = dir(fullfile(cpr_paths().data, 'crlh_far_field.mat'));
fp.dataBytes = d.bytes;
end
