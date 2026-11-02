# Results map

Each row links a result to the scripts that produced it, the main output files, and the time-stamped record
that fixed the analysis before it was run. The first table starts from the manuscript. The tables after it are
organized by the section numbers of the project's internal results document of 2026-10-23, which is not
distributed and serves here only to label rows (its old section number is in brackets where the first
record uses it). Names such as P1, P1′, P2, P3, C, B6, E0–E7 and R10a′ are the names of the analyses in the
time-stamped records; the first table gives their place in the manuscript. Entries such as "old record 7.4" give the
section of the project's earlier results log (before 2026-10-19; not distributed), which describes a result; they are
not time-stamped records.

- **Scripts** are in `code/` unless a path is given. `parts` are the comma-separated stage names the script
  takes as its second argument (see the script docstring). The exact invocations used for the GPU/CPU queues
  are in `run_records/queues/`.
- **Outputs** are relative to the project root `ROOT` (`results/` = summary tables that back the numbers in
  the text; `derived/` = per-domain intermediates). See `docs/DATA.md` for what is deposited.
- **Record** keys refer to `results/prereg_2026-10-17.json` (the first record; round36 = GPU
  plan, round37 = CPU ladder, round39 = scale experiments, round41 = additional experiments, round43 = sensitivities
  added after a result). Analyses added from 2026-10-28 refer to twelve separate records, written after the first record
  was closed, and five corrections: round44 (`prereg_2026-10-28_addendum.json`, `_addendum2.json` …
  `_addendum6.json`; correction `_addendum5_fix1.json`), round45 (`prereg_2026-10-29_addendum7.json`,
  `_addendum8.json`) and round46–48 (`prereg_2026-10-30_addendum9.json` with correction `_addendum9_fix1.json`,
  `_addendum10.json` with corrections `_addendum10_fix1.json` and `_addendum10_fix2.json`, `_addendum11.json`) and round49
  (`prereg_2026-10-31_addendum12.json` with correction `_addendum12_fix1.json`). The date in the file names of the
  round45–49 records is Korea Standard Time. Local copies of these and of the
  round34 and round37 entries are in `registration/`.
- Every timestamp is UTC.

## From the manuscript to the rows

| Manuscript | Scripts | Rows |
|---|---|---|
| Framed chains and the edit distance: theorems (S1 Text) and the general lemma on synthetic sequences (S1 Fig) | `r41_synth.py`, `r41_synth_adv.py` | 1.1, 1.2 |
| Coordinate round trip, certification of 500 pairs and computational cost (Results; S2 Text 1, 6, 7) | `r39_roundtrip.py`; `r38_s0.py` (b0, the first check of S2 Text 1 against a grid of 10⁵ rotations), `r38_static.py`, `r39_sec6.py` (with the share within 0.002 of the best value found, not specified in advance, first computed in the notebook cell of `docs/INTERACTIVE_CHECKS.md` 14); `r41_compute_prep.py`, `r41_compute.py` | 2.4–2.6, 2.2 |
| Mean point independent of the starting frame (Methods) | `r45_mpstart.py` | Round45 |
| TM-align's objective: 12 test pairs, factorial comparison and its sensitivity to d₀, linear-form counterexample, chaining, ranking of 1,674 pairs (Fig 2; S2 Text 8) | `docs/INTERACTIVE_CHECKS.md` 1–2; `pu_tm_factorial_v2.py`, `pu_tm_d0_sensitivity.py` and `docs/INTERACTIVE_CHECKS.md` 8 (violation rates); `framespace_wip/experiments/sippl_linear_counterexample.py` and `sippl_identity_check.py` (its comparison with the identity correspondence); `chaining.py`, `c4_chain_extra.py`; `hu1_rulers.py`, `dfrechet_1674.py`, `hu1_kkl_measures.py` (ranking) | 2.1 |
| Calibration on heating: values, size exponent, detection probability, temperature resolution (Fig 3; S3 Text 4); single-domain resolution of +28 K | `r39_D_stem.py`, `r39_core.py` (e0–e4; the 450 K version without the missing value, `results/r39_E4_nancheck.json`, from the notebook cell in `docs/INTERACTIVE_CHECKS.md` 13); `e0_ht_extra.py` (round29) and `e0_resolution.py` (round32) for the single-domain resolution | 3.6, 3.3 |
| Table 2: replicate and frame audits; trajectory window and the diagnostic of its difference from the earlier audit (added after the result); survival filter | `r44_inv.py`; `r39_core.py` (e2), with the earlier audit of `r37_trunk.py` and the diagnostic in a notebook cell, `docs/INTERACTIVE_CHECKS.md` 10 (`results/r39_E2_aligned_diag.json`); notebook cell in `docs/INTERACTIVE_CHECKS.md` 9 (`results/r37_size_survival.json`) | Round44; 3.6, 3.1; 3.1 |
| Folded-versus-unfolded comparisons (Results; S3 Text 4) | `r39_D_unfold.py`, `r39_e6.py` | 3.7 |
| Comparison across sequences (B6): rank correlation with TM-score and its size control, η₊, far-limit prediction, W₂ (Fig 4A) | `r38_b6.py` (including the certified static distances between medoids used by the far-limit prediction); `r44_b6size_cell26.py`; `r44_w2p1.py` | 4.1, Round44 |
| κ of all cells and of the size-matched cells (Fig 4B, S2 Table); spread share of κ² | `r43_sens.py` (s1), from the per-domain values of rows 3.4 (`r38_dff_D.py`, the whole-trajectory force-field reference), 3.6, 4.1, 4.6, 4.9 and 4.11; `r44_inv.py` (split) | 3.8, Round44 |
| Bond-level comparison: d against aligned coordinates on the same sequence (E7), sources of the compression, contribution of the reduced computation | `r38_e7_30.py`, `r38_e7_100.py`; `r41_dcomp.py` and `docs/INTERACTIVE_CHECKS.md` 5; `r44_redsame.py`; `r38_s0.py` (bias) | 3.4, 3.5, Round44, 2.2 |
| Errors of a temperature-conditioned generator (P2; Fig 5A): reading, secondary criterion and its interval, sensitivity analyses; direction and profiles (E5; Fig 5B); folded counts | scripts of row 4.11 (`r38_p2_D.py`, `r38_p2_aggregate_D.py` and the generation scripts); `r43_sens.py` (s4); `r38_compare.py` (aligned-coordinate run means and per-bond spreads of the cells in Fig 5B), `r39_core.py` (e5); `r45_folded_frac.py` | 4.11, 4.10, Round45 |
| Thermophilic and mesophilic homologs: θ* test (P1), its coverage simulations (Methods) and its diagnostic; replication U2 (P1′; Fig 6) with the folded-fraction control; absolute windows; composition of the pairs | `r38_p1.py` (two-sided coverage in part sim) and notebook cells in `docs/INTERACTIVE_CHECKS.md` 11 (one-sided coverage, condition 2), `r38_p1_diag.py`; `r38_p1n.py`; `r44_p1win.py`; `r44_p1comp2.py` | 4.13, 4.14, Round44 |
| Sequence identity: first test (P3) and its exploratory interval; confirmatory set (C) with C-1, C-2, the network-model sensitivity, the nonlinear control and the shorter-domain controls; bond-direction profiles; pair weights (Fig 7) | `r38_p3.py`, `r38_p3_explore2.py`; `r38_c.py`, `r43_sens.py` (s3), `r38_c_nl.py`; `r45_f1.py`; `r44_pairw.py` | 4.4, 4.5, 3.8, Round44, Round45 |
| Table 3: other correspondences, Gly/Pro, agreement with TM-align on medoids; decomposition of the unshared part | `r44_corr_v2.py`, `r44_corr2.py`; `docs/INTERACTIVE_CHECKS.md` 4 | Round44, 4.5 |
| Alignment accuracy shared by identity and r_dyn; identity range (S3 Table) | `r46_align.py`, `r46_align_fix1.py`; `r48_identity_range.py` | Round46–48 |
| Ties in rank-based estimates (Methods); chain breaks; reproduction by framedchain 0.4.1 | `r47_ties.py`, `r47_combine_fix2.py`; `docs/INTERACTIVE_CHECKS.md` 6, `r49_chainbreak_fix1.py`; `docs/INTERACTIVE_CHECKS.md` 7, `r44_repro.py` | Round46–48, Round49, Round44 |
| Upward bias of the reduced computation and widening of the band, incomplete U-statistic, alignment uncertainty, truncated chains and the diagnostics of their failed check (S2 Text 2–5) | `r38_s0.py` (bias, ustat, unc, trunc), `r38_s0_bias2.py`; notebook cells in `docs/INTERACTIVE_CHECKS.md` 12 (diagnostics) | 2.2, 2.3 |
| Negative type of d and sensitivity of η₊ (S3 Text 2) | `code/negtype/*.js`; commands and printed output in `docs/INTERACTIVE_CHECKS.md` 15 | 1.1 |
| Figures 1–7 and S1 Fig; S1 Table; S2 and S3 Tables | `make_figs.py`; `make_s1.py`; `make_tables.py` | Figures and tables |

## Result 1 — the space and the distance

| Section | Result | Scripts (parts) | Main outputs | Record |
|---|---|---|---|---|
| 1.1 | Metric property, normalized distance D̂, TM-align objective = D_TM at fixed d0 | proofs in the manuscript; `framespace_wip/stage0/rulers.py` (`tm_align`, `tm_cost`), tests `framespace_wip/tests_stage0/test_rulers.py` (README, Quick start) | — | not an analysis (proofs in S1 Text; package tests) |
| 1.1 | D is not of negative type (length-2 structures, SO(3) Fourier coefficient, J D J eigenvalue) | `code/negtype/negtype.js`, `so3fourier.js`, `len2check.js` (Node.js); η₊ for random two-bond chains: `negtype_eta_len2.js` | console output (values in S3 Text 2) | not specified in advance: computations run 2026-10-18 07:54–07:56 UTC, coefficients up to l = 14 recomputed with SciPy on 2026-10-30 (output 14:51:48), η₊ for two-bond chains on 2026-10-31 (output 08:55:09), times from the session log; commands and output in `docs/INTERACTIVE_CHECKS.md` 15 |
| 1.2 [3.6r] | General lemma outside proteins (synthetic examples) | `r41_synth.py`, `r41_synth_adv.py` | `results/r41_synth.json`, `results/r41_synth_adv.json` | round41 g6, 2026-10-21 08:31:54 (`r41_synth.py`, random generator); round41 g6_adv, 21 October 09:25:10, added after the g6 result (`r41_synth_adv.py`, adversarial generators, S1 Fig) |

## Result 2 — existing measures and the computation

| Section | Result | Scripts (parts) | Main outputs | Record |
|---|---|---|---|---|
| 2.1 [2] | Implementation and tests | `framespace_wip/stage0/`, `framespace_wip/tests_stage0/` | — | not an analysis (implementation and its tests) |
| 2.1 [2] | Factorial comparison of TM-score forms (blocks A, B), d0 sensitivity | `pu_tm_factorial.py`, `pu_tm_factorial_v2.py` (block B, args: setup json, output tsv), `pu_tm_d0_sensitivity.py` | `results/pu_tm_factorial_*` (block A), `results/pu_tm_factorial_B_pairs.tsv`, `results/pu_tm_factorial_B_d0sens.tsv`; violation rates of block B (Fig 2A) in `results/pu_tm_factorial_B_triangle.tsv`, written by the notebook cell in `docs/INTERACTIVE_CHECKS.md` 8 | not specified in advance (16 October, before the records began; Drive times of the outputs: block A pairs 15:34–16:32 and rates 16:45:01, block B pairs 16:33–20:04, d₀ sensitivity 20:05–21:53, rates of Fig 2A 20:20:19; the criterion was written in working notes, old record 7.3–7.4) |
| 2.1 [2] | Sippl linear form is not a metric (constructed counterexample) | `framespace_wip/experiments/sippl_linear_counterexample.py`; test `test_linear_form_is_not_a_metric`; comparison of the (a, c) values with the identity correspondence: `sippl_identity_check.py` (added 2026-11-02) | result in the script header; the header's "0.3% 안" (within 0.3%) compares with the identity correspondence at the identity superposition (`tm_cost(P, Q, I, E, z, d0)`), for which the printed gaps are 0.25–0.77%; the comparison in S2 Text 8.3, the identity correspondence with its optimal superposition (equal to within 0.001%), is printed by `sippl_identity_check.py` (output in its header) | not specified in advance (17 October 05:31, before the records began; old record 7.4A); `sippl_identity_check.py`: not specified in advance (2026-11-02) |
| 2.1 [2] | Ranking of existing measures, 1,674 pairs; discrete Fréchet; KKL measures | `hu1_rulers.py`, `dfrechet_1674.py`, `hu1_kkl_measures.py` | `results/hu1_rulers_1674.tsv`, `results/dfrechet_1674.tsv`, `results/hu1_kkl_measures.json`; the Spearman correlation of the discrete Fréchet distance with D_TM/√L̄ (0.566, S2 Text 8.6) is taken from the old results log, because the code of that step was not kept (`dfrechet_1674.py` writes the distances only, and `hu1_kkl_measures.py` does not read them) | not specified in advance (the ranking was run on 16–17 October, before the records began; `results/hu1_kkl_measures.json` was written on 17 October at 05:29 UTC; old record 7.6A) |
| 2.1 [2] | Fold retrieval, 208 structures (not reported in the manuscript) | `pu_fold208.py`, `analysis_fold208.py`, `analysis_fold208_ablation.py` | `results/pu_fold208_pairs.tsv`, `results/pu_fold208_retrieval.json`, `results/pu_fold208_ablation.json` | not specified in advance (16–17 October, before the records began: `pu_fold208_retrieval.json` 2026-10-17 01:33:17, `pu_fold208_ablation.json` last modified 06:14:49, Drive times); not reported in the manuscript (old record 7.7) |
| 2.1 [2] | Chain linking, 150 two-domain chains (Fig 2B, S2 Text 8.4) | `chaining.py`, `c4_chain_extra.py` (v2) | `results/chaining.json`, `results/chaining_v2.json`, `derived/chaining.tsv` | `round17_chaining` (2026-10-18 00:20:07), `round29_review_math` (08:40:18) |
| 2.1 [2] | TM-align comparison (12 pairs) and reversed chains (8 structures) | computed interactively with `stage0.rulers` (no standalone script was kept); cell code in `docs/INTERACTIVE_CHECKS.md` 1–2 | printed output in `docs/INTERACTIVE_CHECKS.md` 1–2 (no output file); tables in old record 7.2, 7.5 | not specified in advance (2026-10-16, before the records began; cells of 15:24–16:41 UTC, session log) |
| 2.2 [3.6d] | First session: correctness of the reduced computation | `r38_s0.py` (pick, b0, e7, trunc, ens, bias, ustat, unc), `r38_s0_bias2.py` | `results/r38_s0_<part>.json`, `results/r38_s0_bias2.json` | round36 §0 (full text in `round36_partial_freeze_7b7c`, 2026-10-19 14:56:51); settings and choice rule of the re-measurement in `r38_s0_bias2.py`: `round36_trackB_code_2026-10-19` (19:43:30; first entered 19:42:45; added after the result) |
| 2.3 [3.6d] | Truncated chains | `r38_s0.py` (trunc); diagnostics of the failed check (S2 Text 5): two scripts written and run by notebook cells, `docs/INTERACTIVE_CHECKS.md` 12 | `results/r38_s0_trunc.json`, `derived/r38/s0/trunc/`; diagnostics in `results/r38_s0_trunc_diag.json` and `results/r38_s0_trunc_exact.json` | round36 §0.1; diagnostics not specified in advance (cells of 19 October, 16:38–16:56 UTC) |
| 2.4 [3.6l] | Coordinate round trip | `r39_roundtrip.py` | `results/r39_roundtrip.json` | `round36_roundtrip_2026-10-20` (stamp 17:23:15; reached Drive 17:34:46, after the output of 17:26:04, as flagged in S1 Table), which also names the carbonyl template and the mean-point projection as descriptive |
| 2.5 [3.6n, 3.6t] | Branch-and-bound certification, static axis 500 pairs | `r38_static.py` (dom, align, bnb, stats), `r39_sec6.py` (time, stats) | `derived/r38/static/bnb.jsonl`, `results/r38_static.json`, `results/r39_sec6.json` | round36 §6 (`round36_sec6_2026-10-20`, 2026-10-20 17:34:46) and §7g (`round36_amend_C_2026-10-20`, 06:23:45); the share within 0.002 of the best value found (499/500) is not specified in advance: notebook cell after the result (output 19:05:44, session log; `docs/INTERACTIVE_CHECKS.md` 14), definition in `round36_sec6` fix1 (19:06:13) and part stats of `r39_sec6.py` (`bnb_posthoc`) |
| 2.6 [3.6w] | Compute table | `r41_compute_prep.py`, `r41_compute.py` | `results/r41_compute.json`; equipment check in `docs/INTERACTIVE_CHECKS.md` 3 | round41 g5 (`g5_design`, 2026-10-22 07:33:57; `g5_impl` 07:41:19) |

## Result 3 — ensemble comparison and the unit

| Section | Result | Scripts (parts) | Main outputs | Record |
|---|---|---|---|---|
| 3.1 [3.1] | Unit invariance: replicate/frame count, drift, trajectory length, size, survival, R0 band | `r37_trunk.py` (R0 band, R1 second half), `ladder_half.py` (drift), `c1_core.py`; the R1 size and survival versions were computed in a notebook cell from the per-domain files of `ladder_temp_single_v2.py` and `c1_core.py` (`docs/INTERACTIVE_CHECKS.md` 9) | `results/r37_trunk.json`, `results/r37_size_survival.json`, `results/ladder_half.json`, `derived/c1_core/` | round37 R0, R1 (2026-10-18 17:08:49); round34 ③ (15:16:42) |
| 3.2 [3.2] | Lower anchor: same-sequence temperature response (folded frames) | `mdcath_5T_extract.py`, `mdcath_5T_folded.py`, `c1_core.py`, `ladder_temp_single.py`, `ladder_temp_single_v2.py`, `ladder_components.py`, `ladder_half.py` | `results/mdcath_5T_setup.json`, `results/ladder_temp_single.json`, `results/ladder_components.json`, `derived/ladder_components_temp.tsv`, `results/ladder_half.json` | rounds 27, 32, 34 (2026-10-18 05:50:23, 12:52:14, 15:16:42; old record 0.2) |
| 3.3 [3.6] | Resolution (E0) | `e0_ht_extra.py` (round29: spread criterion, energy distance; its part 3, the H-T check, is not used in the manuscript), `e0_resolution.py`, `c1_core.py` | `results/e0_ht_extra.json`, `results/e0_resolution.json`, `derived/e0_resolution.tsv` | `round32_review3` E0_resolution (2026-10-18 12:52:14; outcome 13:34:51) and `round29_review_math` E0 (08:40:18; outcome 08:49:09), old record 0.2 |
| 3.4 [3.6d, 3.6f] | D and aligned coordinates on the same sequence (E7) | `r38_s0.py` (e7), `r38_e7_30.py`, `r38_e7_100.py` (pick, run, sum), `r38_dff_D.py` | `results/r38_s0_e7.json`, `results/r38_e7_30.json`, `results/r38_e7_100.json`, `results/r38_dff_D.json` | round36 §3 (full text in `round36_partial_freeze_7b7c`, 2026-10-19 14:56:51); §7g(a), (b) (`round36_amend_C_2026-10-20`, 2026-10-20 06:23:45; script checksums 06:51:03), added after the E7 30 result |
| 3.5 [3.6s] | Mechanism of the compression of D | `r41_dcomp.py` (cpu, gpu parts) | `results/r41_dcomp.json`; where the shifts occur (per-domain `derived/r41/dcomp_gpu/`): `docs/INTERACTIVE_CHECKS.md` 5 | round41 g2 (2026-10-21 08:31:54), added after the E7 result; the counts in end bonds and core are its GPU part, summarized in the cell of `docs/INTERACTIVE_CHECKS.md` 5 (14:45:57, session log) |
| 3.6 [3.6t] | Scale experiments E0′–E4 | `r39_core.py` (e0, e1, e2, e3, e4), `r39_D_stem.py` (run, sum), `r39_sec6.py` | `results/r39_E0.json` … `results/r39_E4.json`, `results/r39_D_stem.json`, `results/r39_sec6.json` | `round39_2026-10-20` (2026-10-20 17:34:46), `round39_amend1_2026-10-21` (2026-10-21 08:27:34); the 450 K check without the missing value (`results/r39_E4_nancheck.json`) is not specified in advance: notebook cell, output 12:37:48 (session log; `docs/INTERACTIVE_CHECKS.md` 13) |
| 3.6 | Diagnostic of the window difference: E2 rerun on the aligned-coordinate values of the same cells (added after the E2 result) | notebook cell in `docs/INTERACTIVE_CHECKS.md` 10 (runs part e2 of `r39_core.py` with the D values replaced by the aligned-coordinate values) | `results/r39_E2_aligned_diag.json` | `round39_posthoc_window_diag_2026-10-21` (2026-10-21 12:43:00 UTC, entered by the same cell before the computation) |
| 3.7 [3.3, 3.6v] | Upper anchor: unfolding (R7) and its position (E6) | `r37_unfold.py`; `r39_D_unfold.py` (run, sum), `r39_e6.py` | `results/r37_unfold.json`, `results/r39_unfold_D.json`, `results/r39_E6.json` | round37 R7 (2026-10-18 17:08:49); round39 E6 in `round39_amend1_2026-10-21` (2026-10-21 08:27:34; analysis script 09:36:39) |
| 3.8 [3.6x] | Length matching of the κ cells (Fig 4B, S2 Table), interval of the P2 secondary criterion, network-model sensitivity of C-2 (added after a result) | `r43_sens.py` | `results/r43_sens.json` | `round43_2026-10-22` s1, s3, s4 (2026-10-22 08:54:20), added after the results |

## Result 4 — sequence differences

| Section | Result | Scripts (parts) | Main outputs | Record |
|---|---|---|---|---|
| 4.1 [3.6i] | Continuous axis (B6), with the certified static distances between the two medoids of each pair used by the far-limit prediction (branch-and-bound in part ens, `static` in each per-pair file) | `r38_b6.py` (ens, stats); check `r38_verify_b.py` | `results/r38_b6.json`, `derived/r38/b6/`, `results/r38_verify_b.json` | round36 §1; rules in `r38_b6.py`, checksum in `round36_trackB_code_2026-10-19` (2026-10-19 19:43:30; first entered 19:42:45) |
| 4.2 [3.6n] | Static axis: medoid distances of the 12,003 H-pool pairs and of the 4,140 of the 6,000 candidate cross-superfamily pairs whose domains had ensembles (128 starts), and branch-and-bound certification of 500 random pairs among them (S2 Text 1); not the distances of the far-limit prediction of 4.1, which `r38_b6.py` certifies itself | `r38_static.py` | `results/r38_static.json` | round36 §1, §7g (`round36_amend_C_2026-10-20`, 2026-10-20 06:23:45; `r38_static.py` checksum 06:51:03) |
| 4.3 [3.6q] | Claim 1: static-axis split and information in the gap share | `r41_gapinfo.py` | `results/r41_gapinfo.json` | round41 g1 (2026-10-21 08:31:54) |
| 4.4 [3.6e, 3.6h] | Primary P3: conservation of dynamics | `r38_p3.py` (pairs, dom, pibar, pair, stats), `r38_p3_sens.py`; exploratory, not specified in advance: `r38_p3_explore.py`, `r38_p3_explore2.py` | `results/r38_p3.json`, `results/r38_p3_sens.json`, `results/r38_p3_explore*.json` | round36 §1b; rules in `r38_p3.py`, checksum in `round36_p3_spec_2026-10-19` (2026-10-19 18:37:38); the exploratory scripts are not specified in advance, their results entered afterwards in `round36_amend_C_2026-10-20` (2026-10-20 06:23:45) |
| 4.5 [3.6j, 3.6x] | Confirmatory experiment C | `r38_c_probe.py`, `r38_c_fetch.py`, `r38_c.py` (pairs, dom, anm, pibar, pair, anmpair, stats), `r38_verify_c.py`; robustness check recorded after the result: `r38_c_nl.py` | `results/r38_c.json`, `results/r38_verify_c.json`, `results/r38_c_nl.json`; decomposition of the unshared part between mean points (descriptive, not specified in advance): `docs/INTERACTIVE_CHECKS.md` 4 | round36 §1c (`round36_amend_C_2026-10-20`, 2026-10-20 06:23:45); `round36_c_nl_2026-10-20` (15:26:31; synthetic check in its amend1, 15:31:34), added after the result; the decomposition is not specified in advance (cells of 2026-10-21 09:17–09:23, session log; values entered afterwards in `round42_withdrawn_2026-10-21`, 09:27:41) |
| 4.6 [3.6o] | R10a′: natural insertions/deletions | pairs selected in a notebook cell from the TM-align and D_TM insertion/deletion events of `p2_gap_consensus.py` (`results/p2_gap_consensus_raw.tsv`) and stored as `r10a_prime` in `results/r37_trunk_pairs.json`; `r38_r10a.py` (ens, stats); added after the result, `r39_r10a_posthoc.py` | `results/r38_r10a.json`, `results/r39_r10a_posthoc.json` | round36 §2 (criteria in `round36_r10a_2026-10-20`, 2026-10-20 16:55:13, stamp 16:43:20); post hoc: `round36_r10a_posthoc_2026-10-20` (19:33:49) |
| 4.7 [3.6o] | R6: insertions/deletions in the same protein (descriptive) | `r39_r6.py` (ens, stats) | `results/r39_r6.json` | round36 §2; `round36_r6_2026-10-20` (2026-10-20 17:34:46) |
| 4.8 [3.4] | Branch cells (same sequence, different source): force field + protocol, protocol only, model error, substitution | `ladder_ff.py`, `ladder_ff_folded.py`, `ladder_ff_single.py`, `ladder_model_error.py`, `ladder_model_error_folded.py`, `ladder_components.py`, `c4_substitution.py`, `c4_chain_extra.py`, `r37_ff.py`, `r38_dff.py`, `r37_protocol.py`, `r37_subst.py` | `results/ladder_ff*_summary.json`, `results/ladder_model_error*_summary.json`, `results/ladder_components.json`, `results/c4_substitution.json`, `results/c4_ladder_row.json`, `results/r37_ff.json`, `results/r38_dff.json`, `results/r37_protocol.json`, `results/r37_subst.json` | rounds 9, 10, 14, 23, 27, 29, 32 (old record 0.2); round37 R3′, R3″, R5′ (2026-10-18 17:08:49) |
| | *Note for 4.8:* the `ladder_ff*` and `ladder_model_error*` scripts write per-domain files (`derived/ladder_ff*/`, `derived/ladder_model_error*/`). The `*_summary.json` files that hold the medians in the table were aggregated in notebook cells that were not kept as scripts; `ladder_components.py` recomputes the component columns from the same per-domain files. | | | |
| 4.9 [3.6u] | R3′ window-matched d_ff | `r41_dff50.py` | `results/r41_dff50.json` | round41 g4 (2026-10-21 08:31:54) |
| 4.10 [3.6m] | E5 direction test | `r39_e5prep.py` (r3p, r3pp), `r39_core.py` (e5) | `results/r39_E5.json` | `round39_2026-10-20` E5 (2026-10-20 17:34:46) |
| 4.11 [3.6g, 3.6x] | Primary P2: ordering of aSAMt 320 K errors, main version D | `r38_prep.py`, `r38_md_ref.py`, `r38_md_ref_p90.py`, `r38_asamt_gen.py` (samenv), `r38_alphaflow_gen.py` (afenv), `r38_compare.py`, `r38_compare_450.py`, `r38_p2_aggregate.py`, `r38_p2_D.py`, `r38_p2_aggregate_D.py`, `r38_p2_D_af.py`, `r38_p2_aggregate_D2.py` | `results/r38_prep.json`, `results/r38_md_ref.json`, `derived/r38/{asamt,alphaflow,compare,compare_D,compare_D2}/`, `results/r38_p2_<TAG>.json` (e.g. `r38_p2_p2_67_D2.json`) | round36 §7c (design and secondary criteria in `round36`, 2026-10-18 18:58:51; full text with the primary reading rule in `round36_partial_freeze_7b7c`, 2026-10-19 14:56:51), §7d (`round36_amend_p2_D_2026-10-19`, 17:07:53); D version of the secondaries in `round36_amend_P2_D2_2026-10-20` (2026-10-20 11:25:00) |
| 4.12 [3.5] | H-T scalar version (θ*) | `mdcath_spread.py`, `mdcath_spread_hot.py`, `ht_theta_meso.py`, `folded_analyses.py` (iii), `ht_own_alignment.py` | `results/ht_shared_spread.tsv`, `results/ht_theta_star.tsv`, `results/ht_theta_meso.json`, `results/folded_analyses.json`, `results/ht_own_alignment.json` | old record 6A; round34 ① (2026-10-18 15:16:42) |
| 4.13 [3.6k] | Primary P1: thermal adaptation θ* | `r38_p1.py` (pairs, mean, align, spread, sim, pre, stats), `r38_verify_p1.py`; one-sided coverage and condition-2 simulations before the decision: two scripts written and run by notebook cells, `docs/INTERACTIVE_CHECKS.md` 11; diagnosis after the result: `r38_p1_diag.py` | `results/r38_p1.json`, `results/r38_p1_sim.json` (two-sided coverage, part sim), `results/r38_p1_sim_onesided.json`, `results/r38_p1_sim_cond2.json`, `results/r38_verify_p1.json`, `results/r38_p1_diag.json` | round36 §7 (estimator simulation and its coverage rule in the freeze checklist of `round36_partial_freeze_7b7c`, 19 October 14:56:51, and in `round36_p1_code_2026-10-20`, 04:58:59, stamp 04:58:52, which also holds the checksum of `r38_p1.py` and so pins part sim); the one-sided simulation applies that rule; the condition-2 simulation (false-positive rate and power of condition 2) is not specified in advance; the results of both were entered in `round36_p1_checklist_2026-10-20` before the decision file was read (one-sided 12:57:39, stamp 12:57:26; condition 2 12:58:12, stamp 12:58:00; flagged in S1 Table); diagnosis: `round36_p1_diag_2026-10-20` (added after the result) |
| 4.14 [3.6p] | Confirmatory P1′ on new pairs | `r38_p1n_fetch.py`, `r38_p1n.py` (mean, align, spread, fold, stats), `r38_verify_p1n.py` | `results/r38_p1n.json`, `results/r38_verify_p1n.json` | round36 §7h (`round36_amend_P1prime_2026-10-20`, 2026-10-20 13:36:03), including the descriptive U1, Q ≥ 0.75 and 0.80, 50–100 ns at Q ≥ 0.75 and the part outside the homolog pool |

## Round44 — checks added after external review (2026-10-28)

Each row except 2.7 and 2.8 (not specified in advance) was specified after the result it qualifies was known and before it was
computed (manuscript, Methods).

| Section | Result | Scripts (parts) | Main outputs | Record |
|---|---|---|---|---|
| 2.7 | framedchain 0.4.1 reproduces the as-run code on synthetic inputs | `r44_repro.py` | printed output in `docs/INTERACTIVE_CHECKS.md` 7 (no output file) | not specified in advance (implementation check; runs 2026-10-28 05:14 and 05:15) |
| 2.8 | No chain breaks (Cα–Cα > 4.2 Å) in the calibration and homolog-study domains | cell code in `docs/INTERACTIVE_CHECKS.md` 6 | printed output there (no output file) | not specified in advance (data check; cells of 2026-10-28 05:15:39–05:22:39, session log); the other sets in Round49 |
| 3.9 | Invariance audit, Table 2 (replicates, frames, split) | `r44_inv.py` (rep, frm, split, sum) | `results/r44_inv.json` | addendum 04:47:45 `r44_inv` |
| 4.1 | Continuous axis controlled for log n̄ | `r44_b6size_cell26.py` | `results/r44_b6size.json` | addendum3 06:04:10 `r44_b6size` |
| 4.5 | Pair-weighted estimates of C-1, C-2 and the H-pool | `r44_pairw.py` | `results/r44_pairw.json` | addendum2 05:56:13 `r44_pairw` |
| 4.5 | Alternative correspondences (TM-align on medoids, sequence alignment) and Gly/Pro | `r44_corr_v2.py` (corr 8 k, stats); first version `r44_corr.py` (stopped) | `results/r44_corr.json` (with `r44_gp`) | addendum 04:47:45 `r44_corr`; addendum2 05:56:13 `r44_corr_v2`, `r44_gp` |
| 4.5 | Agreement against TM-align with tolerances; TM-align on mean points | `r44_corr2.py` (corr 8 k, stats); not specified in advance: `r44_seqid_cell30.py` | `results/r44_corr2.json` | addendum4 06:46:46 `r44_agree`, `r44_tmmean` (including the common pair set with its own seed, the second bootstrap of Table 3) |
| 4.14 | P1′ in absolute windows (50–100, 300–400, 400–500 ns) | `r44_p1win.py` (spread, stats) | `results/r44_p1win.json` | addendum2 05:56:13 `r44_p1win` (the 50–100 ns window at Q ≥ 0.75 is also a pre-specified descriptive item of P1′, row 4.14) |
| 4.1 | Source of the gap between W₂ and Δ²_D (spread) | `r44_w2p1.py` (w2) | `results/r44_w2gap.json` | addendum5 08:53:17 `r44_w2gap` |
| 3.5 | Contribution of the reduced computation to the compression of d on same-sequence cells | `r44_redsame.py` (run, sum) | `results/r44_redsame.json` | addendum6 11:50:14 `r44_redsame` |
| 4.14 | Composition of thermophile and control pairs; U2 adjusted for it | `r44_p1comp2.py` (first attempt: `r44_w2p1.py` p1comp, stopped) | `results/r44_p1comp.json` | addendum5 `r44_p1comp`, corrected by `_addendum5_fix1.json` 08:59:08 |

## Round45 — checks added after external review (2026-10-28 and 10-29)

| Section | Result | Scripts (parts) | Main outputs | Record |
|---|---|---|---|---|
| 4.5 | Sequence-identity effect with bond-direction (tilt) and twist profiles; agreement of network-model predictions with each | `r45_f1.py` (corr 8 k, sanity, stats; queue `q_r45.sh`) | `results/r45_f1.json` | addendum7 2026-10-28 15:39:50 UTC `r45_f1` |
| 2.6 | Mean point independent of the starting frame (six ensembles × nine starts) | `r45_mpstart.py` | `results/r45_mpstart.json` | addendum8 2026-10-29 11:15:02 UTC `r45_mpstart` |
| 4.11 | Share of folded MD 320 K frames and folded aSAMt samples in the P2 domains | `r45_folded_frac.py` | `results/r45_folded_frac.json` | not specified in advance (descriptive count; run on Colab 2026-10-29 11:31–11:37 UTC, output written 11:36:55) |

## Round46–48 — checks added after external review (2026-10-29 and 10-30)

| Section | Result | Scripts (parts) | Main outputs | Record |
|---|---|---|---|---|
| 4.5 | Whether alignment accuracy shared by identity and r_dyn accounts for their association: identity over TM-align's aligned residues on medoids (A) and the precision of the edit-distance correspondence as a covariate, exact and within two bonds (B0, B2) | `r46_align.py` (tm, stats; the stats run stopped at its second check), `r46_align_fix1.py` (average ranks, as corrected) | `derived/r46/align_pairs.jsonl`, `derived/r46/tm_check.json`, `results/r46_align.json`, `results/r46_align_fix1.json` | addendum9 2026-10-29 16:34:05 `r46_align`, corrected by `_addendum9_fix1.json` 17:19:37 |
| Methods | Ties in rank-based estimates whose intervals end within 0.011 of zero: ranks from the original sort against average ranks | `r47_ties.py` (ties; p3, nl, anm × argsort, average; combine), run by `run_records/queues/q_r47.py`; comparison step re-run by `r47_combine_fix2.py` (`r47_combine_fix1.py` stopped with a syntax error) | `derived/r47/<part>_<mode>.json`, `results/r47_ties.json` | addendum10 2026-10-30 06:13:48 `r47_ties`, corrected by `_addendum10_fix1.json` 06:46:19 and `_addendum10_fix2.json` 06:47:24 |
| 4.5 | Whether the association extends below 15% identity: C-1, C-2 below 15% and 10%, expected values under range restriction, residual r_dyn by identity decile, composition-expected identity | `r48_identity_range.py` (exp, stats) | `derived/r48/exp_ident.jsonl`, `results/r48_identity_range.json` | addendum11 2026-10-30 07:30:01 `r48_identity_range` |

## Round49 — chain-break check of further sets (2026-10-31)

| Section | Result | Scripts (parts) | Main outputs | Record |
|---|---|---|---|---|
| Methods | No Cα–Cα bond longer than 4.2 Å in the sets not checked on 2026-10-28: 67 model-error domains, 494 domains of the comparison across sequences, 34 of the natural insertion–deletion pairs, 570 of the certified pairs, 659 and 531 of the two thermophile studies, 37 force-field chains and 666 static chains (2,594 checks: 1,693 distinct mdCATH domains, 37 force-field chains and 666 static chains; longest bond 4.149 Å; the sets not checked are listed in S1 Table) | `r49_chainbreak_fix1.py` (the recorded `r49_chainbreak.py` stopped without output when the Colab VM was replaced) | `results/r49_chainbreak.json`, `derived/r49/items.jsonl` | addendum12 2026-10-31 13:51:43 `r49_chainbreak`; fix1 14:25:18 |

## Figures and tables

The figure and table scripts read stored results only and recompute nothing. On Colab, all eight figures were last
drawn together on 2026-11-02 at 06:07 UTC, the S1 and S3 Tables on 2026-11-02 at 06:07 UTC and the S2 Table on
2026-11-01 at 18:29 UTC (written again with the same content on 2026-11-02 at 06:07 UTC; times of the files on Colab; earlier runs are listed in `docs/PROVENANCE.md`); the scripts
write to `manuscript_figs/` under the project path set at the top of each file.

| Item | Script (function) | Reads |
|---|---|---|
| Fig 1 | `make_figs.py` (fig1) | — (schematic) |
| Fig 2 | `make_figs.py` (fig2) | `results/pu_tm_factorial_B_triangle.tsv` (written by the notebook cell in `docs/INTERACTIVE_CHECKS.md` 8), `results/chaining.json`, `results/chaining_v2.json` |
| Fig 3 | `make_figs.py` (fig3) | `results/r39_E0.json` … `r39_E3.json`, `derived/r39/D/*.json` |
| Fig 4 | `make_figs.py` (fig4) | `derived/r38/b6/*.json`, `results/r43_sens.json`, `derived/r44/split.json` |
| Fig 5 | `make_figs.py` (fig5) | `results/r38_p2_p2_67_D.json`, `results/r39_E5.json` |
| Fig 6 | `make_figs.py` (fig6) | `results/r38_p1n.json`, `results/r44_p1win.json` |
| Fig 7 | `make_figs.py` (fig7) | `results/r38_p3.json`, `r38_p3_explore2.json`, `r38_c.json`, `r38_c_nl.json`, `r44_pairw.json`, `r44_corr.json`, `r44_corr2.json`, `r45_f1.json` |
| S1 Fig | `make_figs.py` (s1fig) | `results/r41_synth_adv.json` |
| S1 Table | `make_s1.py` | `results/_registration_evidence_2026-10-26/entries_timeline.tsv` and the added records `results/prereg_2026-10-*_addendum*.json`; the rows of the predictions and of the classification of analyses are written in the script, transcribed from the records and result files |
| S2 Table | `make_tables.py` (s2) | `results/r43_sens.json` |
| S3 Table | `make_tables.py` (s3) | `results/r46_align_fix1.json`, `results/r48_identity_range.json` |

## Box — standard-metric check (7b)

| Section | Result | Scripts (parts) | Main outputs | Record |
|---|---|---|---|---|
| 5.1 [3.6b] | Reproduction of the AlphaFlow table | `r38_af_bar.py` (base environment), `environment/af_bar_readmeenv.sh` (afenv) | `results/r38_af_bar.txt`, `results/r38_af_bar_readmeenv*.txt` | round36 §7b 0 (predictions in `round36`, 2026-10-18 18:58:51; full text in `round36_partial_freeze_7b7c`, 2026-10-19 14:56:51) |
| 5.2 [3.6c] | Verdict and model ranking | `r38_7b_split.py` (afenv), `r38_7b_rank.py`, `r38_7b_judge.py` | `derived/r38/b7/`, `results/r38_7b_rank.json`, `results/r38_7b_summary.json` | round36 §7b (rule frozen 14:56:51 UTC 2026-10-19) |

## Shared building blocks

These are used by many rows above and are loaded by other scripts with `exec` (keep the `code/` layout):

- `r38_s0.py` — reduced computation, ensemble statistic, sampling helpers (loaded by `r38_b6`, `r38_e7_*`,
  `r38_p2_D`, `r39_core`, `r41_*`).
- `r38_p3.py` — MD readers, folded-frame filter, medoids (loaded by `r38_p1`, `r38_c`, `r38_static`,
  `r38_p3_sens`).
- `r38_compare.py` — model vs MD comparison in aligned coordinates (loaded by `r38_compare_450`).
- `c1_core.py` — ladder CPU core on `derived/mdcath_5T` (loaded by `ladder_temp_single*`, `ladder_components`,
  `e0_resolution`).
- `ladder_model_error.py` — readers and statistics shared by the `ladder_ff*` scripts.
- `framespace_pkg/framespace` (framespace 0.2.0) and `framespace_wip/stage0` — imported via `sys.path`.
