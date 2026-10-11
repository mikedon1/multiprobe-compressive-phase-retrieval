# Multi-Probe Compressive Phase Retrieval: simulation code and data

MATLAB code and data for the simulations in

> M. Don and G. R. Arce, "Compressive Phase Retrieval for Multi-Probe Antenna Measurement Systems"
> (submitted to IEEE Transactions on Antennas and Propagation).

The code recovers antenna radiation patterns from amplitude-only (phaseless) multi-probe
measurements. Random sets of probes are combined with on/off switches and calibrated cables, several
tones are transmitted in each reading, and the pattern is recovered with a sparsity prior and
PR-GAMP. It also contains the comparison with the two-scan phaseless near-field technique
(Section V) and the full-wave example with HFSS probe responses (Section VI).

## Requirements

* MATLAB (developed and tested with R2025b on Windows). **No toolboxes are needed.**
* **GAMPmatlab**, which is *not* included in this repository (it carries no license, so it cannot be
  redistributed). Download it from <https://sourceforge.net/projects/gampmatlab/> (release
  `gampmatlab20171011` or the trunk) and tell the code where it is, in one of three ways:
  set `S.gampPath` in `matlab/cpr_settings.m`, set the `GAMPMATLAB` environment variable, or unpack it into
  `external/GAMPmatlab`. The code uses `code/main` and `code/phase` of that package.

## Quick start

```matlab
>> setup                 % once per MATLAB session, from the repository root
>> check_reproduction    % about 30 s: reproduces one number of the paper (IS, m/n = 0.3: -43.7 dB)
>> run_all('Preset', 'quick')   % everything at reduced size (see below)
```

The quick run writes to `results_quick/` and took about 1.5 hours on a laptop (10 physical cores). `run_all` with no argument uses the paper settings and takes days. `run_all` lists in its comments what each step produces and which paper figure it belongs to.
`run_study('pr_ratio')` runs a single study, and `cpr_load('pr_ratio')` returns its results as a table.

## Settings

Everything you are likely to change is at the top of **`matlab/cpr_settings.m`**:

| Setting | Meaning | `paper` | `quick` |
|---|---|---|---|
| `nw`, `nd` | pattern resolution (theta samples; phi uses `nw/2`) and number of wavelengths | 64, 16 | 32, 4 |
| `restarts` | maximum random restarts, [IS/RS, IC/RC] | [10 30] | [3 6] |
| `nitEM`, `nitGAMP`, `solverSNR` | EM updates, GAMP iterations per update, stopping setting SNR_s | [6 20], [200 1000], [40 48] | smaller |
| `trials` | trials per point (sweeps, circular SNR sweep, near-field, full-wave) | 1, 3, 1, 10 | 1, 1, 1, 2 |
| `sweeps` | the x-axis values of every sweep | the paper's | a few points |
| `nearfield` | SNR and maximum number of two-scan starts | 50 dB, 400 | 50 dB, 10 |
| `fullwave` | number of HFSS frequencies and GAMP iterations in Section VI | 16, 200 | 4, 100 |
| `figures.print` | false: only draw figures in MATLAB; true: also write PDF and PNG | true | true |
| `figures.axes` | `'paper'`: axis limits and ticks as in the paper; `'auto'`: automatic, from the data | `'paper'` | `'auto'` |

`cpr_settings('quick')` or `cpr_settings('paper')` (or `run_all('Preset', 'quick')`) switches the preset for the current session. The other
base parameters (4 active probes, 4 tones per reading, cable lengths U(0, 50) cm, ...) are in
`matlab/core/cpr_cfg.m`. The quick preset exists to check that everything runs and to see the trends fast;
its numbers do **not** match the paper.

## Repository layout

```
matlab/core/              readings, sensing operator, sparse basis, solver, scoring
matlab/studies/           study definitions, loading, conjugate-ambiguity checks
matlab/nearfield/         Section V: spherical-wave expansion, near-field data, two-scan technique
matlab/fullwave/          Section VI: HFSS data, leg calibration, full-wave trials
matlab/figures/           figure scripts
data/                     antenna far field (CRLH leaky-wave antenna, openEMS) and HFSS Touchstone files
```

## Notes on reproducibility

* Every job seeds MATLAB's random generator from its seed, so reruns are exact on the same MATLAB version.
* With the `paper` preset, Sections IV and V use 1000 GAMP iterations per EM update for IC and RC, and
  Section VI uses 200 (`S.fullwave.nitGAMP`), which is what the reported Section VI runs used.
* Everything is checkpointed per job: an interrupted run resumes where it stopped, and
  `run_study(name, 'Workers', N, 'Worker', k)` splits a study over several MATLAB processes.

## Data

* `data/crlh_far_field.mat`: far field of the CRLH leaky-wave antenna, simulated with openEMS
  (T. Liebig et al.) using the parameters of the openEMS CRLH leaky-wave antenna tutorial. Only the fields
  the code uses are kept (E_theta and E_phi at the first 32 frequencies, 1.90 to 3.45 GHz).
* `data/hfss/`: Touchstone files exported from Ansys HFSS-IE (2016.2) for the Section VI example:
  32 rotations of a 20-port model (AUT dipole, 16 ring probes, 3 leg probes) and the reference-dipole solve
  used for the leg calibration.

## License and citation

The code and data in this repository are released under the MIT License (see `LICENSE`). GAMPmatlab is a
separate package with its own terms. If you use this repository, please cite the paper and the archived
release (see `CITATION.cff`).
