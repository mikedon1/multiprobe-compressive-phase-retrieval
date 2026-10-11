function [Xp, info] = load_hfss_mf(files, varargin)
%LOAD_HFSS_MF  Multi-frequency "true pattern" from per-rotation HFSS exports.
%
%   [Xp, info] = load_hfss_mf(globOrFiles, 'AutPort', 1, 'PhiVar', 'AUT_psi')
%
%   Xp   : nProbe x nPhi x nFreq complex AUT->probe transmission, one row per
%          non-AUT port IN FILE ORDER, columns sorted by rotation angle.
%          (load_hfss_pattern.m does the same for the first frequency only.)
%   info : .psi_deg (1 x nPhi), .freq (1 x nFreq, Hz), .portNames, .files
%
% The caller splits ring vs leg rows by the VERIFIED port order
% (cpr_fw_check_ports.m), never by assumption.
p = inputParser;
p.addParameter('AutPort', 1);
p.addParameter('PhiVar', 'AUT_psi');
p.parse(varargin{:});
aut = p.Results.AutPort; phiVar = p.Results.PhiVar;
if ischar(files) || isstring(files)
    d = dir(char(files));
    if isempty(d), error('load_hfss_mf:noFiles', 'No files match %s', files); end
    files = fullfile({d.folder}, {d.name});
end
files = cellstr(files);
nF = numel(files); psi = nan(1, nF); cols = cell(1, nF); f0 = [];
for k = 1:nF
    [S, f, pn, ~, vars] = read_touchstone(files{k});
    N = size(S, 1); pr = setdiff(1:N, aut);
    if isempty(f0), f0 = f(:).'; names = pn(pr);
    elseif numel(f) ~= numel(f0) || any(abs(f(:).' - f0) > 1)
        error('load_hfss_mf:freq', '%s has a different frequency list', files{k});
    end
    cols{k} = permute(S(aut, pr, :), [2 1 3]);            % nProbe x 1 x nFreq
    if isfield(vars, phiVar), psi(k) = sscanf(vars.(phiVar), '%f');
    else, error('load_hfss_mf:noPhi', 'no %s in %s', phiVar, files{k}); end
end
[psi, o] = sort(psi);
if numel(unique(psi)) ~= numel(psi), error('load_hfss_mf:dup', 'duplicate rotation angles'); end
Xp = cat(2, cols{o});
info.psi_deg = psi; info.freq = f0; info.portNames = names; info.files = files(o);
end
