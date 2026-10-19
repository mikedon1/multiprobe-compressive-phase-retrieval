function T = read_aedt_table(file)
%READ_AEDT_TABLE  Read an AEDT "Data Table" CSV export with one row per sweep value of the X variable.
%   T = read_aedt_table(file)
% AEDT writes one column per (quantity, other sweep variables), e.g.
%   "Theta [deg]","re(NearEX) [V/m] - Freq='2.375GHz' Phi='0deg'", ...
% T.x       : X values (first column), as written (deg for angles)
% T.xname   : header of the first column
% T.phi     : sorted unique Phi values (deg) found in the column labels ([] if none)
% T.freq    : sorted unique Freq values (Hz) ([] if none)
% T.q.(name): real array [numel(x), numel(phi), numel(freq)] per quantity, name = label before the
%             unit with re(X) -> re_X, im(X) -> im_X; NaN where a combination is absent.
% T.c.(X)   : complex re_X + 1i*im_X when both exist.
fid = fopen(file, 'r');
if fid < 0, error('read_aedt_table:open', 'cannot open %s', file); end
hdr = fgetl(fid); fclose(fid);
cols = regexp(hdr, '"([^"]*)"', 'tokens'); cols = [cols{:}];
if isempty(cols), cols = strsplit(hdr, ','); end
M = readmatrix(file, 'NumHeaderLines', 1, 'Delimiter', ',');
if size(M, 2) ~= numel(cols)
    error('read_aedt_table:cols', '%s: %d header labels, %d data columns', file, numel(cols), size(M, 2));
end
T.x = M(:, 1); T.xname = cols{1};
nc = numel(cols) - 1;
qn = cell(1, nc); fr = nan(1, nc); ph = nan(1, nc); sc = ones(1, nc);
for j = 1:nc
    c = cols{j + 1};
    u = regexp(c, '\[([^\]]*)\]', 'tokens', 'once');   % AEDT picks a unit PER COLUMN (mV, uv, nV ...)
    if ~isempty(u), sc(j) = unit_scale(u{1}); end
    t = regexp(c, '^\s*([^\s\[]+)', 'tokens', 'once');
    qn{j} = regexprep(regexprep(t{1}, '^(re|im|mag|ang_deg|ang_rad)\((\w+)\)$', '$1_$2'), '\W', '_');
    t = regexp(c, 'Freq=''([-+\d.eE]+)\s*(GHz|MHz|kHz|Hz)''', 'tokens', 'once');
    if ~isempty(t)
        fs = struct('GHz', 1e9, 'MHz', 1e6, 'kHz', 1e3, 'Hz', 1);
        fr(j) = str2double(t{1}) * fs.(t{2});
    end
    t = regexp(c, 'Phi=''([-+\d.eE]+)\s*deg''', 'tokens', 'once');
    if ~isempty(t), ph(j) = str2double(t{1}); end
end
T.freq = unique(fr(~isnan(fr))); T.phi = unique(ph(~isnan(ph)));
nf = max(1, numel(T.freq)); np = max(1, numel(T.phi));
names = unique(qn, 'stable');
T.q = struct();
for k = 1:numel(names)
    T.q.(names{k}) = nan(numel(T.x), np, nf);
end
for j = 1:nc
    ip = 1; if ~isnan(ph(j)), ip = find(T.phi == ph(j)); end
    jf = 1; if ~isnan(fr(j)), jf = find(T.freq == fr(j)); end
    T.q.(qn{j})(:, ip, jf) = M(:, j + 1) * sc(j);
end
T.unitScales = unique(sc);
T.c = struct();
for k = 1:numel(names)
    if startsWith(names{k}, 're_')
        b = names{k}(4:end);
        if isfield(T.q, ['im_' b]), T.c.(b) = T.q.(names{k}) + 1i * T.q.(['im_' b]); end
    end
end
end

function s = unit_scale(u)
% SI scale of an AEDT unit label: prefix (f p n u m k M G) on V, V/m, A, A/m, W; others (deg, ...) 1.
u = strtrim(u); s = 1;
base = {'v', 'v/m', 'v_per_meter', 'a', 'a/m', 'a_per_meter', 'w'};
if any(strcmpi(u, base)), return, end
if numel(u) > 1 && any(strcmpi(u(2:end), base))
    p = u(1);
    switch p
        case 'f', s = 1e-15;  case 'p', s = 1e-12;  case 'n', s = 1e-9;
        case {'u', 'U'}, s = 1e-6;  case 'm', s = 1e-3;  case 'k', s = 1e3;
        case 'M', s = 1e6;  case 'G', s = 1e9;
        otherwise, error('read_aedt_table:unit', 'unknown unit prefix in [%s]', u);
    end
end
end
