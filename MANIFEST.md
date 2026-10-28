# Paper item → study → script → results

Each study is defined in `matlab/studies/cpr_jobs.m` and run with `run_study(<name>)` (all of them: `run_all`). Results are
written to the folder set in `cpr_settings.m` (`resultsDir`, or the `CPR_RESULTS` environment variable). `cpr_paper_figs('IV')` and
`cpr_paper_figs('v9')` write the Section VI panels as `pr2_<name>.pdf` and `pr3_<name>.pdf` into `<results>/figs/paper`; the paper uses the
`pr3_` files. Figure and table numbers are those of the submitted paper and supplement (`S` = supplement).

| Paper item | Study or function | Figure script → file |
|---|---|---|
| Fig. 3(a), wavelengths | `pr_nd` | `cpr_paper_figs('IV')` → `pr3_nd.pdf` |
| Fig. 3(b), tones per reading | `pr_K` | `cpr_paper_figs('IV')` → `pr3_K.pdf` |
| Fig. 3(c), active ports | `pr_P` | `cpr_paper_figs('IV')` → `pr3_P.pdf` |
| Fig. 4 and Fig. S5 (azimuth cut) | `pr_examples` | `cpr_cuts_vsh`, `cpr_cuts_vsh('phi')` → `pr3_cut.pdf`, `pr3_cut_phi.pdf` |
| Fig. 5(a), sampling ratio | `pr_ratio` | `cpr_paper_figs('IV')` → `pr3_ratio.pdf` |
| Fig. 5(b), SNR | `pr_snr` | `cpr_paper_figs('IV')` → `pr3_snr.pdf` |
| Fig. 5(c), cable errors and self-calibration | `pr_calerr`, `pr_calerr30`, `pr_calerrsc`, `pr_calerrics`, `pr_calerrrc` | `cpr_paper_figs('v9')` → `pr3_calerr.pdf` |
| Fig. 6, spherical-wave vs. Fourier basis | `pr_basiscmp32` | `cpr_paper_figs('v9')` → `pr3_basiscmp.pdf` |
| Fig. 8 (azimuth cut) and Fig. S6 (elevation cut), worst of the ten HFSS ring trials | `run_fullwave` | `cpr_fwc_cuts` → `cutphi.pdf`, `cuttheta.pdf` |
| Fig. S2, two-scan technique vs. compressive phase retrieval | `twoscan_run`, `twoscan_continue`, `cpr_sphere_run` for IC, IS, RS, RC (`run_nearfield`) | `cpr_fig_sec5` → `ift_freq.pdf` |
| Fig. S3, installed-probe fraction | `pr_rn` | `cpr_paper_figs('IV')` → `pr3_rn.pdf` |
| Fig. S4, IC with the Fourier basis | `pr_icbad` | `cpr_ic_bad_vsh` → `IC_bad.pdf` |
| Fig. S7, acquisition-time model | formula | `cpr_fig_sec5` → `ift_accel.pdf` |
| Sec. VI-A, switches vs. Gaussian weights | `pr_gauss` | text |
| Sec. VI-B, no sparsity prior | `pr_nosparse` | text |
| Sec. IV and Table S4 text, conjugate candidate, span fractions, identifiability margins | `ident_margins`, `ident_spec`, `round7_checks` | printed (`matlab/diagnostics/`) |
| Table S10, control studies | `pr_iters` (iteration budget), `pr_cables` (cable model), `pr_combiner`, `pr_pairs` (port sets), `pr_tonespan` and `pr_jitter` (tone window, jitter), `pr_trunc` and `pr_aut2` (truncation, larger antennas), `pr_caldrift`, `pr_calgain`, `pr_smooth` (smooth cable and gain errors, dispersion), `pr_scprior` (prior width), `pr_repeat` | supplement table, `cpr_load(<study>)` |
| Table S11, finite detector response | `pr_detbw`; group-delay-aware scores: `detector_delay_scores(1)`, `detector_delay_scores(2)` | `matlab/analysis/` |
| Table S12, rotation-phase margin from the amplitude Jacobian | `jac_rot_geometry` | `matlab/diagnostics/` |
| Tables S3 and S4, sampling geometries, identifiability margins | `diag_geo`, `diag_geo2`, `diag_geo3`, `diag_t1c` | `matlab/diagnostics/` (each script calls `setup` itself) |
| Sec. VII-A, model check of the HFSS ring | `fwc_build_ring` (ring calibration −29.24 dB, half-ring −28.12, free-field probes −26.25, opposite rotation sense −0.42) | printed, stored in `data/fullwave/fwc_ring_hfss.mat` (`Tr.check`) |
| Sec. VII-B, ten trials and the m/n 0.6 and perturbed-calibration runs; Tables S5 and S6 | `run_fullwave` (`fwc_run_trials` with `'Ring','true'`, `'RingCal','ring'`, `'Polish',true`) | `summary_fullwave`, `fwc_numbers`, `fwc_joint_delay` in `matlab/analysis/` |
| Sec. VIII, solve times of the two-scan runs | stored two-scan runs | `ts_times`, `cs_times` in `matlab/analysis/` |

Seeds: every sweep point is one trial with seed 1 unless the study says otherwise; the Section VII trials use seeds 1 to 10 (m/n 0.3) and
1 to 3 (m/n 0.6 and the perturbed-calibration runs). Figures 1, 2 and 7 are drawings and a model screenshot.

## Full-wave data (Section VII)

All files are in `data/fullwave/` (compressed MATLAB files). They are processed from HFSS-IE and method-of-moments data; the HFSS projects
and unpacked exports are not included (available from the author on request).

| File | Content | Used by |
|---|---|---|
| `fwc_validation.mat` | free-space AUT mode coefficients (degree-12 fit to the AUT-alone far field), embedded responses of the 32 arch ports (fits to the near fields of the array alone, by reciprocity), solved transmissions of the 33-port arch model | `fwc_run_trials`, `cpr_fwc_cuts`; rebuilt by `fwc_validate_hfss` (`CPR_HFSS`) |
| `fwc_ring_hfss.mat` | transmissions of the true 65-port ring (16 rotations) and the 64-port ring calibration (`Tr`; `Tr.check` holds the model check) | `fwc_ring_data('hfss')`; rebuilt by `fwc_build_ring` (`CPR_HFSS_RING`) |
| `fwc_validation_mom.mat`, `fwc_ring_mom.mat` | thin-wire method-of-moments data (arch and true ring): the accurate-data check | `fwc_ring_data('momval')`, `fwc_ring_data('mom')`; rebuilt with `mom_make_validation`, `mom_make_ring` |
| `hfss_aut.mat`, `ringalone_fullarc_map.csv` | angular grid of the HFSS far-field export, port map | `mom_make_validation`, `fwc_validate_hfss` |

The refinement of the Section VII estimates is `fwc_polish` (option `Polish` of `fwc_run_trials`; the paper uses it for Section VII only).