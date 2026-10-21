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

Not in this repository: the sampling-geometry and identifiability diagnostics of the supplement (Section S5),
the larger-array restarts started at the true pattern (Section VI-C), and the HFSS build scripts. They are
diagnostic scripts run once; ask the author for them.
