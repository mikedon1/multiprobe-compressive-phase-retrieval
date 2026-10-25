# Paper item → study → script → results

Each study is defined in `matlab/studies/cpr_jobs.m` and run with `run_study(<name>)` (all of them: `run_all`).
The results of the paper were computed in several results folders (set with `CPR_RESULTS`); the last column names
the folder that holds the published numbers. `cpr_paper_figs('IV')` writes the Section VI panels as `pr2_<name>.pdf`;
the paper uses them as `pr3_<name>.pdf`.

| Paper item | Study or function | Figure script → file | Results folder |
|---|---|---|---|
| Fig. 3, IC with the Fourier basis | `pr_icbad` | `cpr_ic_bad_vsh` → `IC_bad.pdf` | results_v9 |
| Fig. 4(a), wavelengths | `pr_nd` | `cpr_paper_figs('IV')` → `pr3_nd.pdf` | results_v8 |
| Fig. 4(b), tones per reading | `pr_K` | `cpr_paper_figs('IV')` → `pr3_K.pdf` | results_v8 |
| Fig. 4(c), active ports | `pr_P` | `cpr_paper_figs('IV')` → `pr3_P.pdf` | results_v8 |
| Fig. 5 and supplement azimuth cut | `pr_examples` | `cpr_cuts_vsh`, `cpr_cuts_vsh('phi')` → `pr3_cut.pdf`, `pr3_cut_phi.pdf` | results_v8 |
| Fig. 6(a), sampling ratio | `pr_ratio` | `cpr_paper_figs('IV')` → `pr3_ratio.pdf` | results_v8 |
| Fig. 6(b), SNR | `pr_snr` | `cpr_paper_figs('IV')` → `pr3_snr.pdf` | results_v9 |
| Fig. 6(c), cable errors and self-calibration | `pr_calerr`, `pr_calerr30`, `pr_calerrsc`, `pr_calerrics`, `pr_calerrrc` | `cpr_paper_figs('v9')` → `pr3_calerr.pdf` | results_v9 |
| Fig. 7, spherical-wave vs. Fourier basis | `pr_basiscmp32` | `cpr_paper_figs('v9')` → `pr3_basiscmp.pdf` | results_v8 |
| Supplement, installed-probe fraction | `pr_rn` | `cpr_paper_figs('IV')` → `pr3_rn.pdf` | results_v9 |
| Section VI-A, switches vs. Gaussian weights | `pr_gauss` | text | results_v9 |
| Section VI-B, no sparsity prior | `pr_nosparse` | text | results_v9 |
| Section II, combiner normalization | `pr_combiner` | text | results_v9 |
| Section VI-C, timing jitter, tone window | `pr_jitter`, `pr_tonespan` | text, supplement table | results_v9 |
| Section VI-C, truncation, array antennas | `pr_trunc`, `pr_aut2` | text, supplement table | results_v9 |
| Supplement, self-calibration variants | `pr_scprior`, `pr_repeat`, `pr_caldrift`, `pr_calgain`, `pr_cables`, `pr_pairs` | supplement table | results_v9 |
| Supplement, solver effort | `pr_iters` | supplement table | results_v9 |
| Section IV-C, conjugate checks | `cpr_conj_check`, `cpr_conj_check_k4` | text | printed |
| Section VII, model check | `fwc_validate_hfss` (or the stored `data/fullwave/fwc_validation.mat`) | text | results_sec7 |
| Section VII, ten trials, perturbed calibration | `fwc_run_trials(1:10)`, then `fwc_run_trials(1:3, 'PerturbAzDeg', 0.5)` etc. (`run_fullwave`) | text | results_sec7 |
| Fig. 9 and supplement elevation cut | worst of the ten trials | `cpr_fwc_cuts` → `cutphi.pdf`, `cuttheta.pdf` | results_sec7 |
| Section VIII, Fig. 10 | `twoscan_run`, `twoscan_continue` (two chains), `cpr_sphere_run` for the four architectures (`run_nearfield`) | `cpr_fig_sec5` → `ift_freq.pdf` | results_v9 |
| Supplement, acquisition-time model | formula | `cpr_fig_sec5` → `ift_accel.pdf` | none |

Seeds: every sweep point is one trial with seed 1 unless the study says otherwise; the Section VII trials use seeds
1 to 10 and the perturbed runs seeds 1 to 3. Figures 1, 2 and 8 are drawings and a model screenshot.

Diagnostics of the supplement (Section S5), in matlab/diagnostics/ (run from MATLAB with that folder on the path; each
script changes to the repository root and calls setup itself):

| Supplement item | Script |
|---|---|
| Identifiability margins of the span, rotation phases and IC cable lengths (Table S4) | ident_margins, ident_spec |
| Counterexample fields, structured rotation phases, orthonormalization factor, compressed-regime local tests (Table S4 text) | ound7_checks(part), parts 1 to 5 |
| Margin against grid sampling, role of the pole samples | grid_margin, pole_check |
| Sampling geometry against recovery (one-sided arch, full ring, polar coverage), Section VII pipeline with ideal probes | diag_geo, diag_geo2, diag_geo3, diag_t1c |

The Section VII polish (wc_polish, option Polish of wc_run_trials) is off by default.

Not in this repository (the archive of the paper, which has a DOI, holds them): the raw HFSS exports and projects
of the full-wave model and the HFSS build scripts. The larger-array restarts started at the true pattern
(Section VI-C) are a diagnostic run once with cpr_prgamp options; ask the author for the driver.

Moment-method data for Section VII (`matlab/fullwave/`): `mom_wire_sparams` (thin-wire solver), `mom_ring_wires` (geometry), `mom_make_validation` (arch data: transmissions, AUT modes, embedded responses, model check), `mom_make_ring` (true full ring of 64 probes, 16 rotations over 180 deg, with the embedded responses of the ring alone; use with `fwc_run_trials(..., 'Ring', 'true', 'RingData', Tm, 'RingCal', 'ring' or 'arch')`). The refinement of the Section VII estimates is `fwc_polish` (option `Polish` of `fwc_run_trials`, default off).

`fwc_ring_calibration` builds the 64-port embedded-response calibration of the full ring alone from the near-field exports of the ring-alone HFSS solve (near ports measured, far ports by the exact 180-degree symmetry); use it as `T.Rring` with `fwc_run_trials(..., 'Ring', 'true', 'RingData', T, 'RingCal', 'ring')`.
