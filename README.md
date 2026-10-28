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

* MATLAB (developed and tested with R2025b on Windows). **No toolbox is required**: everything, including GAMPmatlab, runs on base MATLAB (tested
  with the Signal Processing and Optimization Toolboxes removed from the path). The only effect of a toolbox is in Section VIII: with the Optimization Toolbox the
  two-scan benchmark uses `fminunc` (as in the paper), without it a built-in quasi-Newton method and a warning, so the two-scan numbers can differ slightly.
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
| `trials` | trials per point (sweeps, circular SNR sweep, near-field, full-wave) | 1, 1, 1, 10 | 1, 1, 1, 1 |
| `fullwaveExtras` | Section VII: also m/n 0.6, half-ring calibration, perturbed calibration, MoM arch calibration | true | false |
| `sweeps` | the x-axis values of every sweep | the paper's | a few points |
| `nearfield` | SNR and maximum number of two-scan starts | 50 dB, 400 | 50 dB, 10 |
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
matlab/studies/           study definitions (cpr_jobs.m), loading
matlab/fullwave/          Section VII: HFSS and moment-method data, embedded probe responses by reciprocity, full-wave trials
matlab/nearfield/         Section VIII: spherical-wave expansion, near-field data, two-scan technique
matlab/figures/           figure scripts
matlab/diagnostics/       identifiability and sampling-geometry checks of the supplement
matlab/analysis/          scripts that print the numbers quoted in the paper from saved results
data/                     antenna far field (openEMS) and full-wave data (HFSS)
```

## Section VII (full-wave example)

`run_fullwave` runs the true-ring example of the paper: the 65-port HFSS model (AUT plus 64 probe ports, 16 rotations over 180 deg) with the
embedded-probe calibration of the ring alone and the refinement `fwc_polish` (ten trials at m/n 0.3, three at 0.6, four perturbed-calibration
runs, the half-ring calibration of the supplement), the same on the moment-method ring, the cut figures and `summary_fullwave` (the numbers of
Section VII and Tables S5 and S6). It needs only the files in `data/` (below). One trial takes minutes to about an hour.

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
* `data/fullwave/`, the processed full-wave data (Ansys HFSS-IE 2016.2, 16 frequencies from 2.375 to 3.125 GHz). Every file is derived from
  HFSS exports, with the processing described in the supplement and implemented in `matlab/fullwave/` (`fwc_validate_hfss`, `fwc_build_ring`):
  * `fwc_validation.mat`: free-space mode coefficients of the AUT (a degree-12 spherical-wave fit to the AUT-alone far field, residual
    -91.7 dB or better), embedded receiving responses of the 32 probe ports of the one-sided arch (spherical-wave fits to the near fields of the
    array alone on spheres of 12 and 17 cm, converted by reciprocity), and the solved transmissions of the 33-port arch model.
  * `fwc_ring_hfss.mat`: the solved AUT-to-probe transmissions of the true 65-port ring at 16 rotations (entries of the HFSS S-matrix,
    unchanged) and the 64-port ring calibration (the same fit applied to the near fields of the ring alone), double precision; `Tr.check` holds the model check.
  * `fwc_validation_mom.mat`, `fwc_ring_mom.mat`: the same quantities from the thin-wire method of moments (the accurate-data check).
  * `hfss_aut.mat`: the angular grid of the HFSS far-field export (used by `mom_make_validation`); `ringalone_fullarc_map.csv`: port map.
* The original HFSS projects and exports (CSV and Touchstone files, about 1 GB) are not part of the repository; they are available from the
  author on request. `fwc_validate_hfss` (environment variable `CPR_HFSS`) and `fwc_build_ring` (`CPR_HFSS_RING`) rebuild the processed files from them.
  `run_all` and `run_fullwave` use only the processed files.
## License and citation

The code and data in this repository are released under the BSD 3-Clause Clear License (see `LICENSE`). The license grants no rights to any patents; the methods implemented here may be covered by patents of the author. GAMPmatlab is a
separate package with its own terms. If you use this repository, please cite the paper and the archived
release (see `CITATION.cff`).
