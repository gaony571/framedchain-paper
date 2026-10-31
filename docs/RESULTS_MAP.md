# Results map

Each row links a result in the manuscript to the scripts that produced it, the main output files, and the
pre-registration entry that fixed the analysis before it was run. Section numbers follow the project
results document of 2026-10-23 (`atlas/docs/atlas_results.md`); its old section number is in brackets where
the registration record uses it.

- **Scripts** are in `code/` unless a path is given. `parts` are the comma-separated stage names the script
  takes as its second argument (see the script docstring). The exact invocations used for the GPU/CPU queues
  are in `run_records/queues/`.
- **Outputs** are relative to the project root `ROOT` (`results/` = summary tables that back the numbers in
  the text; `derived/` = per-domain intermediates). See `docs/DATA.md` for what is deposited.
- **Registration** keys refer to `results/prereg_2026-10-17.json` (the registration record; round36 = GPU
  plan, round37 = CPU ladder, round39 = scale experiments, round41 = additional experiments, round43 = post-hoc
  sensitivity). Round44 (2026-10-28) keys refer to six separate records, `results/prereg_2026-10-28_addendum.json`, `_addendum2.json` … `_addendum6.json`, and one correction, `_addendum5_fix1.json`, all written after the first record was closed; local copies of these and of the round34 and round37 entries are in `registration/`.
- Every timestamp is UTC.

## Result 1 — the space and the distance

| Section | Result | Scripts (parts) | Main outputs | Registration |
|---|---|---|---|---|
| 1.1 | Metric property, mean normalisation, TM-align objective = D_TM at fixed d0 | proofs in the manuscript; `framespace_wip/stage0/rulers.py` (`tm_align`, `tm_cost`), tests `framespace_wip/tests_stage0/test_rulers.py` (README, Tests) | — | — (theorem/identity) |
| 1.1 | D is not of negative type (length-2 structures, SO(3) Fourier coefficient, J D J eigenvalue) | `code/negtype/negtype.js`, `so3fourier.js`, `len2check.js` (Node.js) | console output, recorded in `atlas_math_lineage.md` §5 | — (computation) |
| 1.2 [3.6r] | General lemma outside proteins (synthetic examples) | `r41_synth.py`, `r41_synth_adv.py` | `results/r41_synth.json`, `results/r41_synth_adv.json` | round41 g6 |

## Result 2 — existing rulers and the computation

| Section | Result | Scripts (parts) | Main outputs | Registration |
|---|---|---|---|---|
| 2.1 [2] | Implementation and tests | `framespace_wip/stage0/`, `framespace_wip/tests_stage0/` | — | — |
| 2.1 [2] | TM factorial design (blocks A, B), d0 sensitivity | `pu_tm_factorial.py`, `pu_tm_factorial_v2.py` (block B, args: setup json, output tsv), `pu_tm_d0_sensitivity.py` | `results/pu_tm_factorial_*`, `results/pu_tm_factorial_B_*` | earlier rounds (old record 7.3–7.4) |
| 2.1 [2] | Sippl linear form is not a metric (constructed counterexample) | `framespace_wip/experiments/sippl_linear_counterexample.py`; test `test_linear_form_is_not_a_metric` | result in the script header | old record 7.4A |
| 2.1 [2] | Ruler map, 1,674 pairs; discrete Fréchet; KKL measures | `hu1_rulers.py`, `dfrechet_1674.py`, `hu1_kkl_measures.py` | `results/hu1_rulers_1674.tsv`, `results/dfrechet_1674.tsv`, `results/hu1_kkl_measures.json` | old record 7.6A |
| 2.1 [2] | Fold retrieval, 208 structures (supplement) | `pu_fold208.py`, `analysis_fold208.py`, `analysis_fold208_ablation.py` | `results/pu_fold208_pairs.tsv`, `results/pu_fold208_retrieval.json`, `results/pu_fold208_ablation.json` | old record 7.7 |
| 2.1 [2] | Chain linking, 150 two-domain chains (supplement) | `chaining.py`, `c4_chain_extra.py` (v2) | `results/chaining.json`, `results/chaining_v2.json`, `derived/chaining.tsv` | prereg round17, round29 |
| 2.1 [2] | TM-align comparison (12 pairs) and reversed chains (8 structures) | computed interactively with `stage0.rulers` (no standalone script was kept); cell code in `docs/INTERACTIVE_CHECKS.md` 1–2 | printed output in `docs/INTERACTIVE_CHECKS.md` 1–2 (no output file); tables in old record 7.2, 7.5 | — |
| 2.2 [3.6d] | First session: correctness of the reduced computation | `r38_s0.py` (pick, b0, e7, trunc, ens, bias, ustat, unc), `r38_s0_bias2.py` | `results/r38_s0_<part>.json`, `results/r38_s0_bias2.json` | round36 §0 |
| 2.3 [3.6d] | Truncated chains | `r38_s0.py` (trunc) | `results/r38_s0_trunc_*.json` | round36 §0.1 |
| 2.4 [3.6l] | Coordinate round trip | `r39_roundtrip.py` | `results/r39_roundtrip.json` | `round36_roundtrip_2026-10-20` |
| 2.5 [3.6n, 3.6t] | Branch-and-bound certification, static axis 500 pairs | `r38_static.py` (dom, align, bnb, stats), `r39_sec6.py` (time, stats) | `derived/r38/static/bnb.jsonl`, `results/r38_static.json`, `results/r39_sec6.json` | round36 §6, §7g |
| 2.6 [3.6w] | Compute table | `r41_compute_prep.py`, `r41_compute.py` | `results/r41_compute.json`; equipment check in `docs/INTERACTIVE_CHECKS.md` 3 | round41 g5 |

## Result 3 — ensemble comparison and the unit

| Section | Result | Scripts (parts) | Main outputs | Registration |
|---|---|---|---|---|
| 3.1 [3.1] | Unit invariance: replicate/frame count, drift, trajectory length, size, survival, R0 band | `r37_trunk.py` (R0 band, R1), `ladder_half.py` (drift), `c1_core.py` | `results/r37_trunk.json`, `results/ladder_half.json`, `derived/c1_core/` | round37 R0, R1; round34 ③ |
| 3.2 [3.2] | Lower anchor: same-sequence temperature response (folded frames) | `mdcath_5T_extract.py`, `mdcath_5T_folded.py`, `c1_core.py`, `ladder_temp_single.py`, `ladder_temp_single_v2.py`, `ladder_components.py`, `ladder_half.py` | `results/mdcath_5T_setup.json`, `results/ladder_temp_single.json`, `results/ladder_components.json`, `derived/ladder_components_temp.tsv`, `results/ladder_half.json` | rounds 27, 32, 34 (old record 0.2) |
| 3.3 [3.6] | Resolution (E0) | `e0_resolution.py`, `c1_core.py` | `results/e0_resolution.json`, `derived/e0_resolution.tsv` | round32 (old record 0.2, round29) |
| 3.4 [3.6d, 3.6f] | D and aligned coordinates on the same sequence (E7) | `r38_s0.py` (e7), `r38_e7_30.py`, `r38_e7_100.py` (pick, run, sum), `r38_dff_D.py` | `results/r38_s0_e7.json`, `results/r38_e7_30.json`, `results/r38_e7_100.json`, `results/r38_dff_D.json` | round36 §3, §7g(b) |
| 3.5 [3.6s] | Mechanism of the compression of D | `r41_dcomp.py` (cpu, gpu parts) | `results/r41_dcomp.json`; where the shifts occur (per-domain `derived/r41/dcomp_gpu/`): `docs/INTERACTIVE_CHECKS.md` 5 | round41 g2 |
| 3.6 [3.6t] | Scale experiments E0′–E4 | `r39_core.py` (e0, e1, e2, e3, e4), `r39_D_stem.py` (run, sum), `r39_sec6.py` | `results/r39_E0.json` … `results/r39_E4.json`, `results/r39_D_stem.json`, `results/r39_sec6.json` | round39, `round39_amend1_2026-10-21` |
| 3.7 [3.3, 3.6v] | Upper anchor: unfolding (R7) and its position (E6) | `r37_unfold.py`; `r39_D_unfold.py` (run, sum), `r39_e6.py` | `results/r37_unfold.json`, `results/r39_unfold_D.json`, `results/r39_E6.json` | round37 R7; round39 E6 |
| 3.8 [3.6x] | Length matching of the Figure 5 cells (post-hoc sensitivity) | `r43_sens.py` | `results/r43_sens.json` | round43 s1 |

## Result 4 — sequence differences

| Section | Result | Scripts (parts) | Main outputs | Registration |
|---|---|---|---|---|
| 4.1 [3.6i] | Continuous axis (B6) | `r38_b6.py` (ens, stats); check `r38_verify_b.py` | `results/r38_b6.json`, `results/r38_verify_b.json` | round36 §1 |
| 4.2 [3.6n] | Static axis, Figure 4 background | `r38_static.py` | `results/r38_static.json` | round36 §1, §7g |
| 4.3 [3.6q] | Claim 1: static-axis split and information in the gap share | `r41_gapinfo.py` | `results/r41_gapinfo.json` | round41 g1 |
| 4.4 [3.6e, 3.6h] | Primary P3: conservation of dynamics | `r38_p3.py` (pairs, dom, pibar, pair, stats), `r38_p3_sens.py`; exploratory, not registered: `r38_p3_explore.py`, `r38_p3_explore2.py` | `results/r38_p3.json`, `results/r38_p3_sens.json`, `results/r38_p3_explore*.json` | round36 P3 |
| 4.5 [3.6j, 3.6x] | Confirmatory experiment C | `r38_c_probe.py`, `r38_c_fetch.py`, `r38_c.py` (pairs, dom, anm, pibar, pair, anmpair, stats), `r38_verify_c.py`; robustness registered after the result: `r38_c_nl.py` | `results/r38_c.json`, `results/r38_verify_c.json`, `results/r38_c_nl.json`; decomposition of the unshared part between mean points (descriptive, not registered): `docs/INTERACTIVE_CHECKS.md` 4 | round36 §1c; `round36_c_nl_2026-10-20` |
| 4.6 [3.6o] | R10a′: natural insertions/deletions | `r38_r10a.py` (ens, stats); post-hoc `r39_r10a_posthoc.py` | `results/r38_r10a.json`, `results/r39_r10a_posthoc.json` | round36 §2 |
| 4.7 [3.6o] | R6: insertions/deletions in the same protein (descriptive) | `r39_r6.py` (ens, stats) | `results/r39_r6.json` | round36 §2 |
| 4.8 [3.4] | Branch cells (same sequence, different source): force field + protocol, protocol only, model error, substitution | `ladder_ff.py`, `ladder_ff_folded.py`, `ladder_ff_single.py`, `ladder_model_error.py`, `ladder_model_error_folded.py`, `ladder_components.py`, `c4_substitution.py`, `c4_chain_extra.py`, `r37_ff.py`, `r38_dff.py`, `r37_protocol.py`, `r37_subst.py` | `results/ladder_ff*_summary.json`, `results/ladder_model_error*_summary.json`, `results/ladder_components.json`, `results/c4_substitution.json`, `results/c4_ladder_row.json`, `results/r37_ff.json`, `results/r38_dff.json`, `results/r37_protocol.json`, `results/r37_subst.json` | rounds 9, 10, 14, 23, 27, 29, 32 (old record 0.2); round37 R3′, R3″, R5′ |
| | *Note for 4.8:* the `ladder_ff*` and `ladder_model_error*` scripts write per-domain files (`derived/ladder_ff*/`, `derived/ladder_model_error*/`). The `*_summary.json` files that hold the medians in the table were aggregated in notebook cells that were not kept as scripts; `ladder_components.py` recomputes the component columns from the same per-domain files. | | | |
| 4.9 [3.6u] | R3′ window-matched d_ff | `r41_dff50.py` | `results/r41_dff50.json` | round41 g4 |
| 4.10 [3.6m] | E5 direction test | `r39_e5prep.py` (r3p, r3pp), `r39_core.py` (e5) | `results/r39_E5.json` | round39 |
| 4.11 [3.6g, 3.6x] | Primary P2: ordering of aSAMt 320 K errors, main version D | `r38_prep.py`, `r38_md_ref.py`, `r38_md_ref_p90.py`, `r38_asamt_gen.py` (samenv), `r38_alphaflow_gen.py` (afenv), `r38_compare.py`, `r38_compare_450.py`, `r38_p2_aggregate.py`, `r38_p2_D.py`, `r38_p2_aggregate_D.py`, `r38_p2_D_af.py`, `r38_p2_aggregate_D2.py` | `results/r38_prep.json`, `results/r38_md_ref.json`, `derived/r38/{asamt,alphaflow,compare,compare_D,compare_D2}/`, `results/r38_p2_<TAG>.json` (e.g. `r38_p2_p2_67_D2.json`) | round36 §7c, §7d; `round36_amend_p2_D_2026-10-19` |
| 4.12 [3.5] | H-T scalar version (θ*) | `mdcath_spread.py`, `mdcath_spread_hot.py`, `ht_theta_meso.py`, `folded_analyses.py` (iii), `ht_own_alignment.py` | `results/ht_shared_spread.tsv`, `results/ht_theta_star.tsv`, `results/ht_theta_meso.json`, `results/folded_analyses.json`, `results/ht_own_alignment.json` | old record 6A; round34 ① |
| 4.13 [3.6k] | Primary P1: thermal adaptation θ* | `r38_p1.py` (pairs, mean, align, spread, sim, pre, stats), `r38_verify_p1.py`; diagnosis after the result: `r38_p1_diag.py` | `results/r38_p1.json`, `results/r38_verify_p1.json`, `results/r38_p1_diag.json` | round36 §7 |
| 4.14 [3.6p] | Confirmatory P1′ on new pairs | `r38_p1n_fetch.py`, `r38_p1n.py` (mean, align, spread, fold, stats), `r38_verify_p1n.py` | `results/r38_p1n.json`, `results/r38_verify_p1n.json` | round36 §7h |

## Round44 — checks added after external review (2026-10-28)

Each row was specified after the result it qualifies was known and before it was computed (manuscript, Methods).

| Section | Result | Scripts (parts) | Main outputs | Registration |
|---|---|---|---|---|
| 2.7 | framedchain 0.4.1 reproduces the as-run code on synthetic inputs | `r44_repro.py` | printed output in `docs/INTERACTIVE_CHECKS.md` 7 (no output file) | not registered (implementation check) |
| 2.8 | No chain breaks (Cα–Cα > 4.2 Å) in the analysed domains | cell code in `docs/INTERACTIVE_CHECKS.md` 6 | printed output there (no output file) | not registered (data check) |
| 3.9 | Invariance audit, Table 2 (replicates, frames, split) | `r44_inv.py` (rep, frm, split, sum) | `results/r44_inv.json` | addendum 04:47:45 `r44_inv` |
| 4.1 | Continuous axis controlled for log n̄ | `r44_b6size_cell26.py` | `results/r44_b6size.json` | addendum3 06:04:10 `r44_b6size` |
| 4.5 | Pair-weighted estimates of C-1, C-2 and the H-pool | `r44_pairw.py` | `results/r44_pairw.json` | addendum2 05:56:13 `r44_pairw` |
| 4.5 | Alternative correspondences (TM-align on medoids, sequence alignment) and Gly/Pro | `r44_corr_v2.py` (corr 8 k, stats); first version `r44_corr.py` (stopped) | `results/r44_corr.json` (with `r44_gp`) | addendum `r44_corr`; addendum2 `r44_corr_v2`, `r44_gp` |
| 4.5 | Agreement against TM-align with tolerances; TM-align on mean points | `r44_corr2.py` (corr 8 k, stats); not registered: `r44_seqid_cell30.py` | `results/r44_corr2.json` | addendum4 06:46:46 `r44_agree`, `r44_tmmean` |
| 4.14 | P1′ in absolute windows (50–100, 300–400, 400–500 ns) | `r44_p1win.py` (spread, stats) | `results/r44_p1win.json` | addendum2 `r44_p1win` |
| 4.1 | Source of the gap between W₂ and Δ²_D (spread) | `r44_w2p1.py` (w2) | `results/r44_w2gap.json` | addendum5 08:53:17 `r44_w2gap` |
| 3.5 | Contribution of the reduced computation to the compression of d on same-sequence cells | `r44_redsame.py` (run, sum) | `results/r44_redsame.json` | addendum6 11:50:14 `r44_redsame` |
| 4.14 | Composition of thermophile and control pairs; U2 adjusted for it | `r44_p1comp2.py` (first attempt: `r44_w2p1.py` p1comp, stopped) | `results/r44_p1comp.json` | addendum5 `r44_p1comp`, corrected by `_addendum5_fix1.json` 08:59:08 |

## Round45 — checks added after external review (2026-10-28 and 10-29)

| Section | Result | Scripts (parts) | Main outputs | Registration |
|---|---|---|---|---|
| 4.5 | Sequence-identity effect with bond-direction (tilt) and twist profiles; agreement of network-model predictions with each | `r45_f1.py` (corr 8 k, sanity, stats; queue `q_r45.sh`) | `results/r45_f1.json` | addendum7 2026-10-28 15:39:50 UTC `r45_f1` |
| 2.6 | Mean point independent of the starting frame (six ensembles × nine starts) | `r45_mpstart.py` | `results/r45_mpstart.json` | addendum8 2026-10-29 11:15:02 UTC `r45_mpstart` |
| 4.11 | Share of folded MD 320 K frames and folded aSAMt samples in the P2 domains | `r45_folded_frac.py` | `results/r45_folded_frac.json` | not registered (descriptive count) |

## Box — standard-metric check (7b)

| Section | Result | Scripts (parts) | Main outputs | Registration |
|---|---|---|---|---|
| 5.1 [3.6b] | Reproduction of the AlphaFlow table | `r38_af_bar.py` (base environment), `environment/af_bar_readmeenv.sh` (afenv) | `results/r38_af_bar.txt`, `results/r38_af_bar_readmeenv*.txt` | round36 §7b 0 |
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
