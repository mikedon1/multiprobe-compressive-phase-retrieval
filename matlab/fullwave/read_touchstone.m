function [S, freq, portNames, opts, vars] = read_touchstone(filename)
%READ_TOUCHSTONE  Minimal Touchstone (.sNp) reader for AEDT/HFSS exports.
%
%   [S, freq, portNames, opts, vars] = READ_TOUCHSTONE(filename)
%
%   Returns:
%     S         N x N x nFreq complex S-parameter matrices
%     freq      nFreq x 1 frequencies, in the file's own units scaled to Hz
%     portNames N x 1 cellstr of terminal names, parsed from the
%               "! Port[k] = name" comment lines AEDT writes (empty if absent)
%     opts      struct with fields .freqUnit .paramType .format .R
%     vars      struct of design variables from AEDT's "!  name = value"
%               header block (e.g. vars.AUT_phi = '45deg'). Names with a
%               leading '$' (project variables) get the '$' stripped and a
%               'g_' prefix, since '$' is not valid in a MATLAB fieldname.
%
%   Deliberately dependency-free (no RF Toolbox), because the MATLAB
%   Student license may not include sparameters(). The AEDT export format
%   is simple and fixed:
%     - option line: "# GHZ S MA R 50.000000"
%     - one data block per frequency: the frequency, then N^2 value pairs
%       in row-major order, WRAPPED across several physical lines
%     - trailing "! Gamma" / "! Port Impedance" comment blocks
%   so we strip comments and stream numeric tokens rather than parsing
%   line-by-line.
%
%   Supported formats: MA (mag/angle-deg), DB (dB/angle-deg), RI (real/imag).

    txt = fileread(filename);
    lines = regexp(txt, '\r\n|\r|\n', 'split');

    opts = struct('freqUnit','GHZ','paramType','S','format','MA','R',50);
    portNames = {};
    vars = struct();
    dataLines = {};

    for k = 1:numel(lines)
        ln = strtrim(lines{k});
        if isempty(ln)
            continue
        end
        if ln(1) == '!'
            % AEDT writes "! Port[1] = port1_T1"
            tok = regexp(ln, '^!\s*Port\[(\d+)\]\s*=\s*(\S+)', 'tokens', 'once');
            if ~isempty(tok)
                portNames{str2double(tok{1}),1} = tok{2}; %#ok<AGROW>
                continue
            end
            % Design variables, e.g. "!  AUT_phi = 45deg" / "!  $AUT_phi = 0deg"
            vtok = regexp(ln, '^!\s*(\$?[A-Za-z_]\w*)\s*=\s*(\S+)\s*$', 'tokens', 'once');
            if ~isempty(vtok)
                vname = vtok{1};
                if startsWith(vname, '$')
                    vname = ['g_' vname(2:end)];   % '$' invalid in fieldnames
                end
                vars.(vname) = vtok{2};
            end
            continue
        end
        if ln(1) == '#'
            parts = regexp(ln(2:end), '\s+', 'split');
            parts = parts(~cellfun(@isempty, parts));
            if numel(parts) >= 1, opts.freqUnit  = upper(parts{1}); end
            if numel(parts) >= 2, opts.paramType = upper(parts{2}); end
            if numel(parts) >= 3, opts.format    = upper(parts{3}); end
            idxR = find(strcmpi(parts, 'R'), 1);
            if ~isempty(idxR) && numel(parts) >= idxR+1
                opts.R = str2double(parts{idxR+1});
            end
            continue
        end
        % strip any trailing inline comment
        ci = strfind(ln, '!');
        if ~isempty(ci)
            ln = strtrim(ln(1:ci(1)-1));
            if isempty(ln), continue; end
        end
        dataLines{end+1} = ln; %#ok<AGROW>
    end

    nums = sscanf(strjoin(dataLines, ' '), '%f');

    % Infer port count from the filename extension (.sNp) when possible,
    % otherwise from the block size.
    [~, ~, ext] = fileparts(filename);
    tok = regexp(ext, '^\.s(\d+)p$', 'tokens', 'once', 'ignorecase');
    if ~isempty(tok)
        N = str2double(tok{1});
    else
        error('read_touchstone:portCount', ...
              'Cannot infer port count from extension "%s"; expected .sNp', ext);
    end

    perFreq = 1 + 2*N*N;
    if mod(numel(nums), perFreq) ~= 0
        error('read_touchstone:badLength', ...
            ['Numeric token count (%d) is not a multiple of %d ' ...
             '(= 1 freq + 2*%d^2 values). File may be malformed.'], ...
            numel(nums), perFreq, N);
    end
    nFreq = numel(nums) / perFreq;

    switch upper(opts.freqUnit)
        case 'HZ',  scale = 1;
        case 'KHZ', scale = 1e3;
        case 'MHZ', scale = 1e6;
        case 'GHZ', scale = 1e9;
        otherwise,  scale = 1;
    end

    S = zeros(N, N, nFreq);
    freq = zeros(nFreq, 1);
    for f = 1:nFreq
        base = (f-1)*perFreq;
        freq(f) = nums(base+1) * scale;
        vals = nums(base+2 : base+perFreq);
        a = vals(1:2:end);
        b = vals(2:2:end);
        switch upper(opts.format)
            case 'MA'
                z = a .* exp(1i*pi/180 .* b);
            case 'DB'
                z = 10.^(a/20) .* exp(1i*pi/180 .* b);
            case 'RI'
                z = a + 1i*b;
            otherwise
                error('read_touchstone:format', ...
                      'Unsupported format "%s"', opts.format);
        end
        % Touchstone is row-major: row i holds S(i,1..N).
        S(:,:,f) = reshape(z, N, N).';
    end

    if numel(portNames) < N
        portNames{N,1} = '';
    end
end
