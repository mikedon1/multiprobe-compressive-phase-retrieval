# Multi-Probe Compressive Phase Retrieval: simulation code and data

MATLAB code and data for the simulations in

> M. Don and G. R. Arce, "Compressive Phase Retrieval for Multi-Probe Antenna Measurement Systems"
> (submitted to IEEE Transactions on Antennas and Propagation).

The code recovers antenna radiation patterns from amplitude-only (phaseless) multi-probe
measurements. Random sets of probe ports are combined with on/off switches and calibrated cables,
several tones are transmitted in each reading, and the pattern is recovered with PR-GAMP and a
structured sparsity prior on band-limited vector spherical harmonics, with a delay-compensated
cosine basis across frequency. It also contains the full-wave example with HFSS-based embedded probe
responses (Section VII), the comparison with the two-scan phaseless near-field technique (Section VIII),
and the self-calibration of cable lengths.

## Requirements

* MATLAB (developed and tested with R2025b on Windows). No toolboxes are needed. The two-scan benchmark
  uses a built-in quasi-Newton method; with the Optimization Toolbox it can use `fminunc` instead.
* **GAMPmatlab**, which is *not* included in this repository. Download it from <https://sourceforge.net/projects/gampmatlab/> (release
  `gampmatlab20171011` or the trunk) and tell the code where it is, in one of three ways:
  set `S.gampPath` in `matlab/cpr_settings.m`, set the `GAMPMATLAB` environment variable, or unpack it into
  `external/GAMPmatlab`. The code uses `code/main` and `code/phase` of that package.

## Quick start

```matlab
>> setup                 % once per MATLAB session, from the repository root
>> check_reproduction    % a few minutes: operator test, then one number of the paper (IS, m/n = 0.3, -46.5 dB)
>> run_all('Preset', 'quick')   % everything at reduced size (see below)
```

`run_all` with no argument uses the paper settings and takes days on a laptop. `run_study('pr_ratio')` runs a
single study, and `cpr_load('pr_ratio')` returns its results as a table. **`MANIFEST.md` maps every figure and
number of the paper to its study, script and results file.**

## Settings

Everything you are likely to change is at the top of **`matlab/cpr_settings.m`**:

| Setting | Meaning | `paper` | `quick` |
|---|---|---|---|
| `nw`, `nd` | pattern grid (theta samples over 0-360 deg; phi uses `nw/2`) and number of wavelengths | 32, 16 | 32, 8 |
| `angBasis`, `vshL` | angular basis and truncation degree | `'vsh'`, 12 | `'vsh'`, 12 |
| `npol` | polarizations measured and recovered | 2 | 2 |
| `restarts` | maximum initializations, [IS/RS, IC/RC] (some studies set their own) | [10 10] | [3 5] |
| `nitEM`, `nitGAMP`, `solverSNR` | EM updates, GAMP iterations per update, stopping setting SNR_s (target -(SNR_s+2) dB) | [6 6], [200 200], [40 40] | [6 6], [200 300], [40 40] |
| `ratioCircular` | sampling ratio m/n of IC and RC | 0.3 | 0.8 |
| `trials` | trials per point (sweeps, circular SNR sweep, near-field, full-wave) | 1, 1, 1, 10 | 1, 1, 1, 2 |
| `sweeps` | the x-axis values of every sweep | the paper's | a few points |
| `nearfield` | SNR and maximum number of two-scan starts | 50 dB, 400 | 50 dB, 10 |
| `fullwave` | number of HFSS frequencies and GAMP iterations in Section VII | 16, 200 | 8, 200 |
| `resultsDir` | results folder (the `CPR_RESULTS` environment variable overrides it) | `results/` | `results_quick/` |
| `figures.print` | false: only draw figures in MATLAB; true: also write PDF and PNG | true | true |
| `figures.axes` | `'paper'`: axis limits and ticks as in the paper; `'auto'`: from the data | `'paper'` | `'auto'` |

`cpr_settings('quick')` or `cpr_settings('paper')` (or `run_all('Preset', 'quick')`) switches the preset for the
current session. The other base values (P = 4 active ports, K = 4 tones per reading, m/n = 0.3, cable lengths
U(0, 50) cm, 100 dB SNR in the sweeps, ...) are in `matlab/core/cpr_cfg.m`. The quick preset exists to check
that everything runs and to see the trends fast; its numbers do **not** match the paper.

## Repository layout

```
matlab/core/              readings, sensing operator, spherical-wave and frequency bases, solver, self-calibration, scoring
matlab/studies/           study definitions (cpr_jobs.m), loading, conjugate-ambiguity checks
matlab/fullwave/          Section VII: HFSS data, embedded probe responses by reciprocity, full-wave trials
matlab/nearfield/         Section VIII: spherical-wave expansion, near-field data, two-scan technique
matlab/figures/           figure scripts
data/                     antenna far field (openEMS) and full-wave data (HFSS)
```

Some files belong to earlier versions of the paper and are kept for reference, not used by `run_all`: the
leg-probe full-wave pipeline (`cpr_fw_trials`, `fw_cal_reference`, `cpr_fw_cuts`, `data/hfss`), the Fourier-basis
figure scripts `cpr_cuts` and `cpr_ic_bad`, the leg-probe studies (`pr_legs`, `pr_legsdft`, `pr_legsmatch32`), and
the complex-measurement studies (`cs_*`), whose results are not in the paper.

## Notes on reproducibility

* Every job seeds MATLAB's random generator from its seed. Reruns reproduce the results on the same
  MATLAB version and platform; numerical libraries and multithreading can change the last digits elsewhere.
* Everything is checkpointed per job: an interrupted run resumes where it stopped, and
  `run_study(name, 'Workers', N, 'Worker', k)` splits a study over several MATLAB processes.
* Each checkpoint stores the settings and data version it was made with. Resuming with different
  settings is refused; use a new results folder instead. Failed jobs are run again on the next call.
* `check_reproduction` raises an error if the operator test or the reference value fails.
* The self-calibration step defaults to the one used for the paper (solve the damped normal equations, then
  remove the mean). `c.selfCalProjectAfter = false` selects the constrained Gauss-Newton step, which changes the
  NMSE by at most 0.15 dB at the points checked.

## Data

* `data/crlh_far_field.mat`: far field of the CRLH leaky-wave antenna, simulated with openEMS
  (T. Liebig et al.) using the parameters of the openEMS CRLH leaky-wave antenna tutorial. Only the fields
  the code uses are kept (E_theta and E_phi at the first 32 frequencies, 1.90 to 3.45 GHz).
* `data/fullwave/fwc_validation.mat`: the processed HFSS data of Section VII (Ansys HFSS-IE 2016.2, 16 frequencies
  from 2.375 to 3.125 GHz): the free-space mode coefficients of the AUT (fit to the AUT-alone solve), the embedded
  receiving response of each of the 32 probe ports to each spherical mode (by reciprocity from the array-alone
  solve), and the solved AUT-to-probe transmissions of the 33-port model at 32 rotations. `run_fullwave` uses it
  directly. `data/fullwave/ringalone_fullarc_map.csv` is the port map.
* The raw HFSS exports (rotation Touchstone files, AUT-alone far field, array-alone near fields on two spheres) are
  needed only to rebuild that file with `fwc_validate_hfss`, which reads them from `data/hfss_fullarc` or from the
  folder named by the `CPR_HFSS` environment variable.
* `data/hfss/`: Touchstone files of the earlier 20-port leg-probe model, not used by the current paper.

## License and citation

The code and data in this repository are released under the BSD 3-Clause Clear License (see `LICENSE`). The license grants no rights to any patents; the methods implemented here may be covered by patents of the author. GAMPmatlab is a
separate package with its own terms. If you use this repository, please cite the paper and the archived
release (see `CITATION.cff`).
