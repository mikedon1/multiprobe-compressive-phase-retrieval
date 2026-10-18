function run_study(study, varargin)
%RUN_STUDY  Run all jobs of one study (one figure's worth of simulations) and checkpoint them.
%
%   run_study('pr_ratio')                       run every job of the study
%   run_study('pr_snr', 'Workers', 4, 'Worker', 2)   run only share 2 of 4 (start several MATLAB
%                                               processes, one per share, to use more cores)
%
% Studies are listed in studies/cpr_jobs.m; sweep grids, trials, restarts and the pattern size come
% from cpr_settings.m. Results go to <results>/<study>_w<Worker>.mat after every job, so an interrupted
% run resumes where it stopped, and a job that is finished in ANY checkpoint of the study is
% never run twice. Collect results with cpr_load(study). Each checkpoint stores the settings it was
% made with (cpr_fingerprint); resuming with different settings is refused, so stale results are never
% mixed with new ones. Failed jobs are recorded but are run again on the next call.
p = inputParser;
p.addParameter('Workers', 1);
p.addParameter('Worker', 1);
p.parse(varargin{:});
NW = p.Results.Workers; WK = p.Results.Worker;
cpr_setup();
S = cpr_settings();
maxNumCompThreads(S.threads);
outDir = cpr_paths().results;
J = cpr_jobs(study);
[~, ord] = sort([J.cost], 'descend');                % greedy longest-first split of the shares
ld = zeros(1, NW); owner = zeros(1, numel(J));
for k = ord, [~, w] = min(ld); owner(k) = w; ld(w) = ld(w) + J(k).cost; end
mine = find(owner == WK);
f = fullfile(outDir, sprintf('%s_w%d.mat', study, WK));
fp = cpr_fingerprint();
if exist(f, 'file'), R0 = load(f); res = R0.res; checkfp(R0, fp, f); else, res = struct([]); end
keyOf = @(j) sprintf('%s|%s|%d', j.study, j.label, j.seed);
doneKeys = okKeys(res);
dd = dir(fullfile(outDir, sprintf('%s_w*.mat', study)));
for i = 1:numel(dd)
    if strcmp(dd(i).name, sprintf('%s_w%d.mat', study, WK)), continue; end
    So = load(fullfile(dd(i).folder, dd(i).name));
    checkfp(So, fp, fullfile(dd(i).folder, dd(i).name));
    doneKeys = [doneKeys, okKeys(So.res)]; %#ok<AGROW>
end
fprintf('%s: share %d/%d, %d jobs (preset %s)\n', study, WK, NW, numel(mine), S.preset);
t0 = tic;
for k = mine
    j = J(k); key = keyOf(j);
    if any(strcmp(doneKeys, key)), continue; end
    c = cpr_cfg(j.over{:});
    try
        r = cpr_run(c, j.seed);
    catch ME
        fprintf('FAILED %s: %s\n', key, ME.message);
        r = struct('cfg', c, 'seed', j.seed, 'time', NaN, ...
            'score', struct('raw', NaN, 'global', NaN, 'perlambda', NaN, 'each', NaN(1, c.nd)), ...
            'info', struct(), 'R', NaN, 'm', NaN, 'n', NaN, 'Xh', single([]));
    end
    e = struct('key', key, 'study', j.study, 'label', j.label, 'seed', j.seed, 'cfg', c, ...
        'raw', r.score.raw, 'global', r.score.global, 'perlambda', r.score.perlambda, ...
        'each', r.score.each, 'time', r.time, 'R', r.R, 'm', r.m, 'info', r.info, 'Xh', []);
    if j.saveXh, e.Xh = r.Xh; end
    if isempty(res), res = e;
    else
        old = find(strcmp({res.key}, key), 1);           % a failed earlier attempt is replaced
        if isempty(old), res(end+1) = e; else, res(old) = e; end %#ok<AGROW>
    end
    save(f, 'res', 'fp');
    fprintf('%s -> raw %7.2f | global %7.2f | per-lambda %7.2f dB (%5.0f s, total %.0f s)\n', ...
        key, e.raw, e.global, e.perlambda, e.time, toc(t0));
end
fprintf('%s: share %d done in %.0f s\n', study, WK, toc(t0));
end

function k = okKeys(res)
% keys of finished jobs; a failed job (NaN score) does not count as finished
k = {};
if isempty(res), return, end
good = arrayfun(@(e) ~(isnan(e.raw) && isnan(e.global)), res);
k = {res(good).key};
end

function checkfp(R, fp, file)
if ~isfield(R, 'fp')
    error('run_study:fingerprint', ['%s has no settings fingerprint (made by an older version). ' ...
        'Use a new results folder (cpr_settings resultsDir or the CPR_RESULTS environment variable).'], file);
end
if ~isequal(R.fp, fp)
    error('run_study:fingerprint', ['%s was made with different settings or data. Use a new results ' ...
        'folder (cpr_settings resultsDir or the CPR_RESULTS environment variable).'], file);
end
end
