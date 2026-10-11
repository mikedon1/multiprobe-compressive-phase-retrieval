% SETUP  Run once per MATLAB session from the repository root:   >> setup
% Adds the code to the path and checks that GAMPmatlab is found (see matlab/core/cpr_setup.m).
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'matlab'), fullfile(root, 'matlab', 'core'));
cpr_setup();
fprintf('Ready. Preset: %s. Try  run_all(''Preset'', ''quick'')  or  run_study(''smoke'').\n', cpr_settings().preset);
clear root
