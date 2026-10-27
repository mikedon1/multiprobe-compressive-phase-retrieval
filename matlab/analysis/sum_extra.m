function sum_extra()
%SUM_EXTRA  Finished jobs of the extra studies (results_extra): label, global NMSE, time.
cpr_setup(); maxNumCompThreads(1);
for f = dir(fullfile(cpr_paths().results, 'pr_*_w*.mat')).'
    S = load(fullfile(f.folder, f.name)); fprintf('== %s: %d jobs\n', f.name, numel(S.res));
    for e = S.res(:).', fprintf('   %-24s global %7.2f  perlambda %7.2f  (%4.0f s)\n', e.label, e.global, e.perlambda, e.time); end
end
end