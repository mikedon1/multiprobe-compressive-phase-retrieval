function T = cpr_load(study)
%CPR_LOAD  Load every worker file of a study into one table.
%   T = cpr_load('cs_map')  -> table with label, seed, meas, arch, P, K, ratio, nw, nd,
%   snrdB, nlegs, rnRS, rnRC, raw, global, perlambda, time, R, m
outDir = cpr_paths().results;
fs = dir(fullfile(outDir, sprintf('%s_w*.mat', study)));
rows = {};
for f = fs.'
    S = load(fullfile(f.folder, f.name));
    for e = S.res(:).'
        c = e.cfg;
        wt = 'switch'; if isfield(c, 'weights'), wt = c.weights; end      % older checkpoints lack the field
        ce = 0; if isfield(c, 'calErrDeg'), ce = c.calErrDeg; end
        rows(end+1, :) = {e.label, e.seed, c.meas, c.arch, c.P, c.K, c.ratio, c.nw, c.nd, ...
            c.snrdB, c.nlegs, c.rnRS, c.rnRC, e.raw, e.global, e.perlambda, e.time, e.R, e.m, wt, ce}; %#ok<AGROW>
end, end
T = cell2table(rows, 'VariableNames', {'label','seed','meas','arch','P','K','ratio','nw','nd', ...
    'snrdB','nlegs','rnRS','rnRC','raw','global','perlambda','time','R','m','weights','calErrDeg'});
% a failed attempt (NaN) recorded in one worker file and rerun successfully in another (e.g. after a change of
% the number of workers): keep only the successful record
key = strcat(T.label, '|', string(T.seed));
bad = isnan(T.raw) & ismember(key, key(~isnan(T.raw)));
T(bad, :) = [];
end
