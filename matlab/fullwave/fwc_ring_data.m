function Tr = fwc_ring_data(kind)
%FWC_RING_DATA  Section VII data of the TRUE full ring, from the files released with the code.
%   Tr = fwc_ring_data('hfss')   65-port HFSS ring: Xring, psi, freq and the ring calibration Rring (data/fullwave/fwc_ring_hfss.mat;
%                                built from the raw HFSS exports by fwc_build_ring, which are not in the repository)
%   Tr = fwc_ring_data('mom')    the same from the thin-wire method of moments (data/fullwave/fwc_ring_mom.mat, mom_make_ring); its
%                                arch validation data (use as 'Validation') come from fwc_ring_data('momval')
if nargin < 1, kind = 'hfss'; end
d = fullfile(cpr_paths().data, 'fullwave');
switch kind
    case 'hfss'
        f = fullfile(d, 'fwc_ring_hfss.mat');
        if ~isfile(f), Tr = fwc_build_ring(f); else, S = load(f); Tr = S.Tr; end
    case 'mom',    S = load(fullfile(d, 'fwc_ring_mom.mat')); Tr = S.Tm;
    case 'momval', S = load(fullfile(d, 'fwc_validation_mom.mat')); Tr = S.V;
    otherwise, error('fwc_ring_data: unknown kind %s', kind);
end
end
