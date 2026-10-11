function P = cpr_paths()
%CPR_PATHS  Locations used by the code: P.root (repository), P.data, P.results.
%   The results folder is <repository>/results (results_quick for the quick preset) unless cpr_settings.m sets resultsDir.
here = fileparts(mfilename('fullpath'));          % <root>/matlab/core
P.root = fileparts(fileparts(here));
P.data = fullfile(P.root, 'data');
S = cpr_settings();
sub = 'results'; if strcmp(S.preset, 'quick'), sub = 'results_quick'; end   % quick runs never mix with paper runs
if isempty(S.resultsDir), P.results = fullfile(P.root, sub); else, P.results = S.resultsDir; end
if ~exist(P.results, 'dir'), mkdir(P.results); end
end
