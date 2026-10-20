function cpr_harvest(resultsDir, studies, outFile)
%CPR_HARVEST  Print every job of the given studies in one table per study, for filling the paper's \TBD values.
%   cpr_harvest('...\results_v9')                       all studies found in the folder
%   cpr_harvest('...\results_v9', {'pr_calerrsc'})      selected studies
%   cpr_harvest(dir, studies, 'harvest_v9.txt')         also write the text to a file
% Per label: number of seeds, NMSE in dB (one global phase for K > 1, one phase per wavelength for K = 1; mean
% over seeds in dB, then min and max when there are several), actual m and R, restarts used (numTry), GAMP runs
% (EM updates + 1, summed over the restarts, when the run logged them), and when present the self-calibration
% length error (rms phase at the band centre, deg: true -> after, common length removed), the gain error (dB),
% and the realized SNR of the combined signals. Reads the worker files directly (cpr_load keeps fewer fields).
% (2026-10-20, round 7)
if nargin < 2 || isempty(studies)
    fs = dir(fullfile(resultsDir, '*_w*.mat'));
    studies = unique(regexprep({fs.name}, '_w\d+\.mat$', ''));
end
if ischar(studies), studies = {studies}; end
txt = {};
for s = studies(:).'
    fs = dir(fullfile(resultsDir, sprintf('%s_w*.mat', s{1})));
    if isempty(fs), txt{end+1} = sprintf('== %s: no files', s{1}); continue, end %#ok<AGROW>
    res = [];
    for f = fs.'
        S = load(fullfile(f.folder, f.name), 'res');
        r = S.res(:);
        r = rmfield(r, intersect(fieldnames(r), {'Xh'}));            % patterns are not needed here
        if isempty(res), res = r; else, res = [res; reshape(orderfields(r, res), [], 1)]; end %#ok<AGROW>
    end
    labels = unique({res.label});
    txt{end+1} = sprintf('== %s (%d jobs, %d labels)', s{1}, numel(res), numel(labels)); %#ok<AGROW>
    txt{end+1} = sprintf('%-28s %3s %8s %8s %8s %6s %5s %5s %6s  %s', 'label', 'n', 'NMSE', 'min', 'max', ...
        'm', 'R', 'try', 'runs', 'extra'); %#ok<AGROW>
    for L = labels
        e = res(strcmp({res.label}, L{1}));
        y = arrayfun(@metric, e); ok = ~isnan(y);
        if any(ok), ym = 10*log10(mean(10.^(y(ok)/10))); else, ym = NaN; end
        nt = arrayfun(@(x) getf(x.info, 'numTry'), e);
        nr = arrayfun(@(x) getf(x.info, 'nRunsTotal'), e);
        ex = '';
        sc = arrayfun(@(x) isfield(x.info, 'selfCal'), e);
        if any(sc)
            t = arrayfun(@(x) getf(getf(x.info, 'selfCal'), 'rmsTrueDeg'), e(sc));
            a = arrayfun(@(x) getf(getf(x.info, 'selfCal'), 'rmsErrDeg'), e(sc));
            ex = sprintf('%s len %.2f->%.2f deg', ex, mean(t, 'omitnan'), mean(a, 'omitnan'));
            gE = arrayfun(@(x) getf(getf(x.info, 'selfCal'), 'rmsGainErrDb'), e(sc));
            if any(~isnan(gE))
                gT = arrayfun(@(x) getf(getf(x.info, 'selfCal'), 'rmsGainTrueDb'), e(sc));
                ex = sprintf('%s gain %.2f->%.2f dB', ex, mean(gT, 'omitnan'), mean(gE, 'omitnan'));
            end
        end
        sn = arrayfun(@(x) getf(x.info, 'snrCombinedDb'), e);
        if any(~isnan(sn)), ex = sprintf('%s SNRc %.1f dB', ex, mean(sn, 'omitnan')); end
        if sum(~ok) > 0, ex = sprintf('%s FAILED %d', ex, sum(~ok)); end
        txt{end+1} = sprintf('%-28s %3d %8.2f %8.2f %8.2f %6d %5d %5.1f %6.1f %s', L{1}, numel(e), ym, ...
            min(y), max(y), e(1).m, e(1).R, mean(nt, 'omitnan'), mean(nr, 'omitnan'), ex); %#ok<AGROW>
    end
end
txt = strjoin(txt, newline);
disp(txt);
if nargin >= 3 && ~isempty(outFile)
    fid = fopen(outFile, 'w'); fprintf(fid, '%s\n', txt); fclose(fid);
end
end

function y = metric(e)
% one phase per wavelength for K = 1 (the phases between wavelengths are not measured), else one global phase
y = e.global;
if isfield(e, 'cfg') && isfield(e.cfg, 'K') && e.cfg.K == 1, y = e.perlambda; end
end

function v = getf(s, f)
v = NaN;
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end
