function [T, c] = fwc_load_trials(pattern)
%FWC_LOAD_TRIALS  Section VII trials saved by fwc_run_trials: all files in <results> whose name matches the regular expression, concatenated
%and sorted by seed. c is the configuration of the first file.
%   [T, c] = fwc_load_trials('^fwc_trials_true_hfss_ring\.mat$')
d = dir(fullfile(cpr_paths().results, 'fwc_trials*.mat')); d = d(~cellfun(@isempty, regexp({d.name}, pattern)));
if isempty(d), error('fwc_load_trials: no file in %s matches %s (run run_fullwave first)', cpr_paths().results, pattern); end
T = []; c = [];
for i = 1:numel(d), S = load(fullfile(d(i).folder, d(i).name)); T = [T, S.T]; if isempty(c), c = S.c; end, end %#ok<AGROW>
[~, o] = sort([T.seed]); T = T(o);
end
