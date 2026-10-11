function cpr_setup()
%CPR_SETUP  Put this repository and GAMPmatlab on the MATLAB path (once per session).
%
% GAMPmatlab is NOT part of this repository. Download it from
%   https://sourceforge.net/projects/gampmatlab/   (release gampmatlab20171011 or the trunk)
% and point to it in ONE of these ways:
%   1. set S.gampPath in cpr_settings.m,
%   2. set the environment variable GAMPMATLAB, or
%   3. unpack it into <repository>/external/GAMPmatlab.
% The folder may be the unpacked project root or its 'code' folder; the code needs the
% 'main' and 'phase' subfolders of 'code' (gampEst, GampOpt, CAwgnEstimIn, SparseScaEstim,
% ncCAwgnEstimOut, ...).
persistent done
if ~isempty(done), return, end
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(genpath(fullfile(root, 'matlab')));
S = cpr_settings();
g = S.gampPath;
if isempty(g), g = getenv('GAMPMATLAB'); end
if isempty(g), g = fullfile(root, 'external', 'GAMPmatlab'); end
codeDir = '';
for c = {g, fullfile(g, 'code'), fullfile(g, 'trunk', 'code'), fullfile(g, 'latest', 'trunk', 'code')}
    if exist(fullfile(c{1}, 'main', 'gampEst.m'), 'file'), codeDir = c{1}; break, end
end
if isempty(codeDir)
    error('cpr_setup:gamp', ['GAMPmatlab not found (looked in "%s"). Download it from ' ...
        'https://sourceforge.net/projects/gampmatlab/ and set gampPath in cpr_settings.m, ' ...
        'or set the GAMPMATLAB environment variable.'], g);
end
addpath(fullfile(codeDir, 'main'), fullfile(codeDir, 'phase'));
if ~exist('ncCAwgnEstimOut', 'file')
    error('cpr_setup:gamp', 'GAMPmatlab found at "%s" but its code/phase folder lacks ncCAwgnEstimOut.m.', codeDir);
end
done = true;
end
