function E = swe_cprh(Et, Ep, phi)
%SWE_CPRH  The openEMS nf2ff E_cprh component (the scalar pattern used throughout the paper):
%   E_cprh = (E_theta + j E_phi) exp(j phi) / sqrt(2)
%(matches CRLH_LeakyWaveAnt_min.mat's E_cprh to -319 dB; cpr_test_swe.m). phi broadcasts.
E = (Et + 1j*Ep) .* exp(1j*phi) / sqrt(2);
end
