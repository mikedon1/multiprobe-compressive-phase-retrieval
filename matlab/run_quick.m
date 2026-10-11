function run_quick(varargin)
%RUN_QUICK  Run everything at reduced size: the 'quick' preset of cpr_settings.m (smaller pattern,
%few wavelengths, few restarts, few sweep points, two full-wave trials). It takes minutes to about
%an hour and writes to results_quick/, so it never mixes with paper-preset results. The numbers
%will NOT match the paper; use it to check that the code runs and to see the trends quickly.
%   run_quick                      all parts
%   run_quick('Skip', {'checks'})  leave parts out (see run_all)
cpr_settings('quick');
run_all(varargin{:});
end