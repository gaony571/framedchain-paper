# Script reference (English)

The docstrings and comments inside the scripts are in Korean, the working language of the project (except the figure
and table scripts `make_figs.py`, `make_tables.py` and `make_s1.py`, which are in English). The
scripts are kept byte-for-byte as run (`PROVENANCE.md` lists their checksums and which of them were recorded before the run), so they are not
translated in place. This page describes every script in English: what it computes, its stages, its main
inputs and outputs, and the record entry it implements. Notation follows the scripts, with the manuscript's terms in
parentheses: an *edge* (bond) is a residue-to-residue step; *D* (d) is the alignment distance of the space (computed
by the reduced computation: banded dynamic programming with mean-point starts, an upper bound); Δ² = 2C − B_P − B_Q
(B_P′ for B_Q) is the energy-distance statistic between two ensembles (C = cross mean distance, B = floors, the mean
distances between different runs within each side, one value per pair of runs); Δ²/√n̄ divides by the square root of
the mean number of edges; Δμ²₂ = (C₂ − (B₂_P + B₂_Q)/2)/edges is the mean-shift component; a *run* (replicate) is one
independent simulation; *folded frames* have Q ≥ 0.75 of the native contacts of the
first 320 K frame unless stated otherwise; *H-pool* is the homolog pool of the manuscript, all 14,189 same-superfamily
pairs of mdCATH domains with at most 150 residues (392 superfamilies).

Usage is `python code/<script>.py ROOT [stages] [workers]` unless stated; stages are resumable.

## Data preparation

| Script | What it does | Main outputs |
|---|---|---|
| `mdcath_5T_extract.py` | Reads each mdCATH h5 once and stores Cα, C, O coordinates (every `frame_step` frames), Q of every frame, Rg, DSSP and the sequence for the 336 domains × 5 temperatures × 5 replicates in `results/mdcath_5T_setup.json` | `derived/mdcath_5T/<dom>.npz` |
| `mdcath_5T_folded.py` | Per-edge frame spread, Cα RMSF² and windowed (±3) RMSF² of the folded frames for each temperature and replicate | `derived/mdcath_5T_folded/<dom>.npz` |
| `mdcath_spread.py`, `mdcath_spread_hot.py` | Per-edge spread of mdCATH domains, frames aligned to the mean point of each temperature. Arguments: temperatures (default 320,348) and an output-directory suffix; the 379/413 K spreads were written with `mdcath_spread_hot.py` to `derived/mdcath_spread_hot` | `derived/mdcath_spread<suffix>/<dom>.npz` |
| `p2_gap_consensus.py` | For the H-pool homolog pairs, lists the insertion/deletion events of the TM-align alignment and of the D_TM alignment (fixed d0 = 4 Å, with the TM-align pose among the starts). The natural insertion–deletion pairs (R10a′) were selected from this table in a notebook cell and stored in `results/r37_trunk_pairs.json` | `results/p2_gap_consensus_raw.tsv` |
| `r38_prep.py` | Test-set preparation for P2: sequences of the 90 aSAMt test domains, sequence identity to the ATLAS training set, AlphaFlow input CSV, ColabFold MSAs | `results/r38_prep.json`, `derived/r38/alphaflow_input.csv`, `derived/r38/msa/` |
| `r38_md_ref.py`, `r38_md_ref_p90.py` | MD reference for P2 domains: mean point of 320 K folded frames, native contacts, atom indices, same-domain temperature curves. `_p90` is a copy that runs only the three domains first excluded for an overlap of PDB IDs (`1i2tA00`, `1nijA01`, `4qbuA03`), completing the 90-domain set, and writes only their per-domain files | `derived/r38/md/<dom>.{npz,json}`, `results/r38_md_ref.json` (`r38_md_ref.py` only) |
| `r38_c_probe.py` | Tests partial HTTP reads of mdCATH h5 files and lists the new-superfamily domains to fetch for experiment C | `results/r38_c_probe.json`, `derived/r38/c_new_sf_domains.json` |
| `r38_c_fetch.py`, `r38_p1n_fetch.py` | Fetch only the needed frames (first frame, second halves; for P1′ also the 50–100 ns window) of new domains from the mdCATH Hugging Face dataset | `raw/mdcath_c/<dom>.npz`, `raw/mdcath_p1n/<dom>.npz` |

## Model ensembles

| Script | What it does | Main outputs |
|---|---|---|
| `r38_asamt_gen.py` | Runs aSAMt (`sam2` commit `ce076ed`, `generate_ensemble.py` with defaults: 250 samples, 100 steps, minimization) per domain and temperature in samenv; seeds fixed per (domain, temperature) by CRC32 because the tool has no seed argument | `derived/r38/asamt/<dom>.<T>.{traj.dcd,top.pdb,time.txt}` |
| `r38_alphaflow_gen.py` | Runs distilled AlphaFlow-MD (commit `0408d7c`, `predict.py` as in its README: 250 samples, `--noisy_first --no_diffusion`) in afenv; seeds fixed per domain | `derived/r38/alphaflow/<dom>.pdb` |

## Comparisons of TM-score forms and existing measures (result 2.1)

| Script | What it does | Main outputs |
|---|---|---|
| `pu_tm_factorial.py`, `pu_tm_factorial_v2.py` | Factorial comparison of TM-score forms: S_TM of the space's optimizer at fixed and length-dependent d0, short/long/mean normalization, number of starts, with and without the TM-align pose, against TM-align (`_v2` takes the setup and output file names as arguments; block B). The triangle-inequality violation rates of block B (Fig 2A, `results/pu_tm_factorial_B_triangle.tsv`) were computed from its pair table in a notebook cell (`docs/INTERACTIVE_CHECKS.md` 8) | `results/pu_tm_factorial_pairs.tsv`, `results/pu_tm_factorial_B_pairs.tsv` |
| `pu_tm_d0_sensitivity.py` | Same pairs with d0 = 2 and 6 Å | `results/pu_tm_factorial_B_d0sens.tsv` |
| `hu1_rulers.py` | Ranking of existing measures on 1,674 pairs: TM-align (both normalizations, RMSD), GDT-TS, lDDT (Cα), S_TM at d0 = 4 Å, D at ℓ = 4, 8, 16, ∞, ESA | `results/hu1_rulers_1674.tsv` |
| `dfrechet_1674.py` | Discrete Fréchet distance of Cα chains minimized over rigid motions (fragment starts + TM-align pose, refined) | `results/dfrechet_1674.tsv` |
| `hu1_kkl_measures.py` | Kolodny–Koehl–Levitt (2005) SI, SAS, MI from the TM-align alignments; rank correlations with the other measures | `results/hu1_kkl_measures.json` |
| `pu_fold208.py` | All pairs of the 208-structure block: TM-align and S_TM at fixed and mean-length d0 | `results/pu_fold208_pairs.tsv` |
| `analysis_fold208.py`, `analysis_fold208_ablation.py` | Topology mean average precision and class 1-NN accuracy, bootstrap difference against TM-align; ablation over normalizations | `results/pu_fold208_retrieval.json`, `results/pu_fold208_ablation.json` |
| `chaining.py` | Chain linking on 150 two-domain chains from CathDomall v4.4 (downloads the PDB files): triangle violations and single/average-link merging for 1 − TM_short, 1 − TM_long, D_TM/√L̄ and D_TM at fixed d0 | `results/chaining.json`, `derived/chaining.tsv` |
| `c4_chain_extra.py` | Round29: substitution ladder row from `derived/c4_sub`, and chain linking v2 with the proven fixed-d0 normalization and linear forms | `results/c4_ladder_row.json`, `results/chaining_v2.json` |
| `framespace_wip/experiments/sippl_linear_counterexample.py` | Constructed counterexample showing that the linear (Sippl) form is not a metric. Its header's "0.3% 안" (within 0.3% of the identity correspondence) refers to the identity correspondence at the identity superposition, for which the printed gaps are 0.25–0.77% (`sippl_identity_check.py`) | result in the file header |
| `sippl_identity_check.py` | Added on 2026-11-02 (not specified in advance) to check the comparison in S2 Text 8.3: for the same chains (seed 0) and the same `stage0.rulers` functions, the optimized (a, c) value against the identity correspondence at the identity superposition and at its optimal superposition (weighted least squares from the identity, 200 rounds); run from the repository root, output in the file header | standard output |

## Ladder, temperature unit and early checks (results 3.1–3.3, 4.8, 4.12)

| Script | What it does | Main outputs |
|---|---|---|
| `c1_core.py` | CPU core on `derived/mdcath_5T`: equilibrium diagnostics, E0, temperature ladder in aligned coordinates (frames of all temperatures aligned to the 320 K folded mean point), exact W2 with POT | `derived/c1_core/<dom>.json` |
| `ladder_temp_single.py`, `ladder_temp_single_v2.py` | Temperature rows with a single-replicate floor and absolute scale (round27); `_v2` adds the floor at each temperature | `results/ladder_temp_single.json`, `derived/ladder_temp_single[_v2]/` |
| `ladder_half.py` | Round34 ③ drift control: first vs second half of trajectories at the same temperature treated as two conditions (expected 0 at equilibrium) | `results/ladder_half.json` |
| `ladder_components.py` | Round32: mean-shift component Δμ²/edge (bilinear floor subtracted) and log total-spread ratio for the temperature, model-error and force-field rows | `results/ladder_components.json`, `derived/ladder_components_temp.tsv` |
| `e0_ht_extra.py` | Round29 E0_and_HT: (1) per-bond spread criterion of E0 on all frames (`derived/mdcath_spread/`): median over bonds of log(mean spread at 348 / 320 K), five replicates each, against the largest of the 20 within-temperature 2:3 replicate splits, with the pass rate and a bootstrap interval; (2) share of domains whose energy distance Δ² at 348 K is positive, from the per-domain files of `ladder_temp_single_v2.py`; (3) H-T check in bacteria with a temperature threshold, which executes part of `c3_saturation_ht_variants.py` (an earlier exploratory script not in this repository; not used in the manuscript), so without that file the script stops at part (3) and writes no output; `results/e0_ht_extra.json` is the output of the Colab run. Written to Google Drive by a notebook cell on 2026-10-18 at 08:44:12 UTC and never modified; restored byte for byte on 2026-11-02. A later cell (08:47 UTC) added a post hoc permutation version, `E0_spread_posthoc_perm`, to the same output | `results/e0_ht_extra.json` |
| `e0_resolution.py` | Round32 E0: permutation test (5 of 10 replicates, 252 splits) of 320 vs 348 K on the same folded frames with three measures: frame energy distance, Cα energy distance, windowed RMSF² log ratio | `results/e0_resolution.json` |
| `ladder_ff.py`, `ladder_ff_folded.py`, `ladder_ff_single.py` | Force-field + protocol row: MFM reference MD ("SM", ff99SB-disp) vs BioEmu CATH MD (amber99sb-ildn) on the same CATH domains; all frames, folded frames, and single-run floors | `derived/ladder_ff[_folded/_single]/<dom>.json` |
| `ladder_model_error.py`, `ladder_model_error_folded.py` | Model-error row: ATLAS test chains (82) vs released AlphaFlow/ESMFlow samples (template and sequence-only); all frames and folded frames. Also provides readers and statistics used by the `ladder_ff*` scripts | `derived/ladder_model_error[_folded]/<chain>.json` |
| `c4_substitution.py` | C4 substitution axis (round14): MegaSim variants vs wild type (ff99SB-disp trajectories only), with Tsuboyama ΔΔG | `derived/c4_sub/<WT>.json`, `results/c4_substitution.json` |
| `ht_theta_meso.py` | Round34 ①: θ* for thermophile–mesophile pairs recomputed, mesophile–mesophile negative control (random-direction and symmetric versions), folded-frame version | `results/ht_theta_meso.json` |
| `folded_analyses.py` | Folded-frame analyses (2026-10-17): (i) Q check, (ii) clean C2 test, (iii) H-T comparison without extrapolation (thermophile 348 K vs mesophile 320 K), (iv) C3 scalar control | `results/folded_analyses.json` |
| `ht_own_alignment.py` | Repeats (iii) with shared edges defined by the space's own D_TM alignment instead of TM-align, and reports the overlap of the two edge sets | `results/ht_own_alignment.json` |

## Round37 CPU ladder (results 3.1, 3.7, 4.8)

| Script | What it does | Main outputs |
|---|---|---|
| `r37_trunk.py` | R0 band (320 K replicate splits) and R1 second-half version of the temperature cells; its whole-trajectory reference values (`full_for_ratio`) are the medians printed by the size and survival audit, which ran in a notebook cell (`docs/INTERACTIVE_CHECKS.md` 9, `results/r37_size_survival.json`) | `results/r37_trunk.json` |
| `r37_ff.py` | R3′: same protein, ATLAS (CHARMM36m, 300 K) vs mdCATH 320 K (CHARMM22*) for 44 chains; temperature-equivalent from the same domain's curve | `results/r37_ff.json` |
| `r37_protocol.py` | R3″: protocol only, BioEmu `ONE_cath1` vs `MSR_cath2` (both amber99sb-ildn) on 16 shared domains | `results/r37_protocol.json` |
| `r37_subst.py` | R5′: substitution Δ² (MegaSim) placed on the same domain's mdCATH temperature curve for 5 systems | `results/r37_subst.json` |
| `r37_unfold.py` | R7 folded vs unfolded: Majewski fast folders (350 K) and mdCATH 413/450 K, in aligned coordinates and d_fix | `results/r37_unfold.json` |

## Round36 GPU plan (results 2.2–2.5, 3.4, 4.1–4.7, 4.11, 4.13–4.14, box 5)

| Script | What it does | Main outputs |
|---|---|---|
| `r38_s0.py` | First session (round36 §0): `pick` domains and pairs with fixed seeds; `b0` branch-and-bound check of the 128-start optimum; `e7` Δ², w̄ and ratio by D and by aligned coordinates on the same samples; `trunc` truncated chains (predicted: the unshared (gap) part of full-versus-truncated comparisons exceeds that of full-versus-full comparisons by 0.75 Σ s² of the removed bonds); `ens` Δ²_D of the 30 pairs; `bias` upper-bound bias of the reduced computation against branch-and-bound by TM bin; `ustat` incomplete U-statistic (2,000 frame pairs) vs the full 100 × 100 mean; `unc` spread of the matched and gap shares across near-optimal alignments (S ≥ (1 − 0.005) S_max). Its first part is loaded by many later scripts | `results/r38_s0_<part>.json`, `derived/r38/s0/` |
| `r38_s0_bias2.py` | Follow-up of §0.2: re-measures the upper-bound bias in the TM 0.2–0.3 bin with more starts and wider bands; decision rule fixed before running | `results/r38_s0_bias2.json` |
| `r38_e7_30.py`, `r38_e7_100.py` | E7: same-sequence temperature cells by D and aligned coordinates on 30 domains (pre-specified; rule frozen 2026-10-19 14:56:51 UTC) and on 100 domains (7g, `round36_amend_C`, 2026-10-20 06:23 UTC, added after the 30-domain result) | `results/r38_e7_30.json`, `results/r38_e7_100.json` |
| `r38_p3.py` | Primary P3 (§1b): conservation of dynamics on H-pool pairs; partial Spearman of r_dyn and sequence identity given the shared-core distance. Stages `pairs, dom, pibar (GPU), pair, stats`. Provides MD readers, folded-frame filters and medoids used by other scripts | `results/r38_p3.json`, `derived/r38/p3/` |
| `r38_p3_sens.py` | Four of the five P3 sensitivity versions (medoid and crystal controls, Q 0.7/0.8, replicate-split reliability correction, log spread ratio) | `results/r38_p3_sens.json` |
| `r38_p3_explore.py`, `r38_p3_explore2.py` | Exploratory verdicts not specified in advance: where the power of P3 is lost; calibration simulation of the pre-specified bootstrap | `results/r38_p3_explore*.json` |
| `r38_b6.py` | B6 continuous axis on 281 trunk pairs: Δ²_D, matched and gap shares, D_matched, fixed-correspondence versions (π̄, TM-align), η₊ of J D J, exact W2, static medoid D | `results/r38_b6.json`, `derived/r38/b6/` |
| `r38_static.py` | Static axis ("Figure 4 background" in the script's working figure numbering, not a figure of the paper): medoid D for H-pool and 6,000 cross-superfamily pairs; branch-and-bound certification of 500 random pairs | `results/r38_static.json`, `derived/r38/static/` |
| `r38_c.py` | Confirmatory experiment C (§1c): P3 repeated on independent superfamilies (C-1) and beyond the dynamics predicted by structure, using ANM edge spreads (C-2) | `results/r38_c.json` |
| `r38_c_nl.py` | Robustness check recorded after the C verdict: flexible (gradient-boosting, cross-fitted) partial correlation, residual shape, general dependence | `results/r38_c_nl.json` |
| `r38_r10a.py` | R10a′ (§2): natural insertion/deletion pairs; pre-specified tests on the R0 band, gap excess vs 0.75 × inserted residues, matched-share residual vs the homolog regression (TOST) | `results/r38_r10a.json` |
| `r38_compare.py` | Common comparison code for P2 in aligned coordinates: all run pairs, Δ², ratio, Δμ², log spread ratio for MD temperature cells, aSAMt and AlphaFlow cells and the MD null | `derived/r38/compare/<dom>.{json,npz}` |
| `r38_compare_450.py` | Adds the aSAMt 450 K cells that were generated after the 67-domain comparison ran (same seeds and call order, so the MD frames are identical) | `derived/r38/compare/<dom>.json` |
| `r38_dff.py`, `r38_dff_D.py` | R3′ force-field difference d_ff with the common code, in aligned coordinates and by D | `results/r38_dff.json`, `results/r38_dff_D.json` |
| `r38_p2_aggregate.py` | P2 aggregation (§7c, frozen; md5 written to the checklist before the run): interval of the aSAMt 320 K error on the same-domain temperature scale, censoring rules, superfamily bootstrap; secondary analyses 1–4 | `results/r38_p2_<TAG>.json` |
| `r38_p2_D.py`, `r38_p2_aggregate_D.py` | P2 by D (amendment `round36_amend_p2_D_2026-10-19`); the aggregate wrapper runs the frozen aggregation on the D inputs without editing it | `derived/r38/compare_D/`, `results/r38_p2_<TAG>.json` |
| `r38_p2_D_af.py`, `r38_p2_aggregate_D2.py` | AlphaFlow cell by D and the aggregation of secondary analyses 2–3 by D | `derived/r38/compare_D2/`, `results/r38_p2_<TAG>.json` |
| `r38_p1.py` | Primary P1 (§7): thermal adaptation θ* on thermophile–mesophile pairs with the symmetric mesophile–mesophile version | `results/r38_p1.json` |
| `r38_p1_diag.py` | Diagnosis after the P1 verdict (does not change it): heterogeneity and power, window comparison, an uncensored paired statistic | `results/r38_p1_diag.json` |
| `r38_p1n.py` | Confirmatory P1′ (§7h) on new pairs: stability vs flexibility (folded fraction by temperature, U2 at Q 0.7/0.8) | `results/r38_p1n.json` |
| `r38_verify_b.py`, `r38_verify_p1.py`, `r38_verify_c.py`, `r38_verify_p1n.py` | Implementation checks of B4/B6, P1, C and P1′ on synthetic data with known answers (and a small real case), without looking at verdict values | `results/r38_verify_*.json` |
| `r38_af_bar.py` | 7b correctness: runs the original AlphaFlow analysis code (commit `0408d7c`) on the four released sample sets and compares with the published table | `results/r38_af_bar.txt` |
| `r38_7b_split.py` | 7b verdict data: exchangeable 10 ns block halves × 20 splits; standard metrics (original code, afenv) and Δ² for the true candidate and 7 controls (shrink, inflate, trim, scramble, crystal + noise) | `derived/r38/b7/` |
| `r38_7b_judge.py` | 7b verdict under the rule frozen at 2026-10-19 14:56:51 UTC: fraction of targets for which each metric calls a control closer than the truth, with target bootstrap | `results/r38_7b_summary.{json,txt}` |
| `r38_7b_rank.py` | Descriptive model ranking: standard-metric medians vs Δ² medians for the four sample sets, Kendall τ | `results/r38_7b_rank.json` |

## Round39 scale experiments (results 2.4, 2.5, 3.6, 3.7, 4.7, 4.10)

| Script | What it does | Main outputs |
|---|---|---|
| `r39_D_stem.py` | Same-sequence cells by D on all 336 trunk domains, with a run-pair table used for the subset statistics | `derived/r39/D/`, `results/r39_D_stem.json` |
| `r39_core.py` | Round39 CPU analyses: `check` (synthetic checks), `pr` (participation ratio), `e0` regime boundary and component ratio, `e1` size exponent α (IRLS), `e2` windows, `e3` power curves, `e4` back-conversion, `e5` direction test | `results/r39_E0.json` … `results/r39_E5.json` |
| `r39_e5prep.py` | Inputs for E5: run means and per-edge spreads of the R3′ and R3″ cells in aligned coordinates | `derived/r39/e5/` |
| `r39_D_unfold.py` | Folded vs unfolded by D on mdCATH 413/450 K with three Q thresholds | `results/r39_unfold_D.json` |
| `r39_e6.py` | E6 position rule (fixed before seeing `r39_unfold_D`): an upper-anchor band is drawn only if the positions of the comparison cells of Fig 4B of the paper ("Figure 5" in the script's working numbering) agree across six bands and under an alternative normalization | `results/r39_E6.json` |
| `r39_r6.py` | R6 (descriptive): three natural insertion/deletion events within the same protein (ATLAS) | `results/r39_r6.json` |
| `r39_r10a_posthoc.py` | Post-hoc analysis after the R10a′ verdict: is the matched-share residual due to the insertions or to other differences between test and reference pairs | `results/r39_r10a_posthoc.json` |
| `r39_sec6.py` | §6 minimal version: timing table of the heuristic, branch-and-bound, TM-align (re-implementation, tmtools) and US-align; certification rate from the 500 static pairs; since its fix1 (round36_sec6, 2026-10-20 19:06 UTC, after the static result), also the post hoc share of pairs within 0.002 of the best value found (`bnb_posthoc`, not specified in advance) | `results/r39_sec6.json` |
| `r39_roundtrip.py` | Coordinate round trip: atoms ↔ coordinates of the space (Cα exact; C, O from a peptide-plane template) and projection error of the mean point | `results/r39_roundtrip.json` |

## Round41 and round43 (results 1.2, 2.6, 3.5, 3.8, 4.3, 4.9)

| Script | What it does | Main outputs |
|---|---|---|
| `r41_gapinfo.py` | Claim 1: static-axis split into matched and gap parts and the information in the gap share (partial Spearman, ΔAUC), and the r_dyn check | `results/r41_gapinfo.json` |
| `r41_dcomp.py` | Mechanism of the compression of D: four levels of distance (aligned coordinates, d_fix, identity correspondence with per-edge gap choice, D) on the same frame pairs as E7 100 | `results/r41_dcomp.json` |
| `r41_dff50.py` | R3′ d_ff with both sides restricted to 50–100 ns | `results/r41_dff50.json` |
| `r41_compute_prep.py`, `r41_compute.py` | Compute table on one machine: TM-align, US-align, Foldseek, D with 1 and 128 starts, branch-and-bound, all-vs-all and ensemble cells, with warm-up and repeat timing | `results/r41_compute.json` |
| `r41_synth.py`, `r41_synth_adv.py` | General lemma outside proteins: triangle inequality on real sequences (ℓ2, ERP), planar curves modulo SO(2) (certified), and negative controls (broken Lipschitz constant, DTW); adversarial triples in the `_adv` version | `results/r41_synth.json`, `results/r41_synth_adv.json` |
| `r43_sens.py` | Round43 post-hoc sensitivities (verdicts unchanged): P2 secondary 1 interval, length matching of the cells of Fig 4B and S2 Table ("Figure 5" in the script's working numbering), ANM/GNM variants for C-2 | `results/r43_sens.json` |

## Round44 and round45 — checks added after external review (2026-10-28 and 10-29; results 2.6, 2.7, 3.5, 3.9, 4.1, 4.5, 4.11, 4.14)

Each analysis was written into a new time-stamped record before it was run (`registration/prereg_2026-10-28_addendum*.json`
and `prereg_2026-10-29_addendum*.json`), except the three checks not specified in advance (`r44_repro.py`, `r44_seqid_cell30.py`,
`r45_folded_frac.py`); the results it qualifies were already known. The first record (`results/prereg_2026-10-17.json`) was not changed.

| Script | What it does | Main outputs |
|---|---|---|
| `r44_inv.py` | Invariance audit of the calibration statistic (Table 2): 3 vs 5 replicates per side from the stored per-pair tables (rep), 50 vs 100 frames per replicate recomputed with D (frm), shares of the split ratio (split), summary (sum) | `results/r44_inv.json` |
| `r44_corr.py` | Confirmatory set C with shared bonds from global sequence alignments and from TM-align on medoids instead of the edit-distance correspondence π̄; superseded after 1,369 pairs by `r44_corr_v2.py` (same formulas) | `derived/r44/corr/pairs.jsonl` (partial, not used) |
| `r44_corr_v2.py` | The same computation in shards (corr K k → stats), plus π̄ without bonds next to Gly/Pro and Gly/Pro covariates (r44_gp); structural control and the network-model profile correlation r_ANM recomputed on each set of shared bonds (the per-domain network-model predictions are reused) | `derived/r44/corr/pairs_v2_s*.jsonl`, `results/r44_corr.json` |
| `r44_pairw.py` | Pair-weighted versions of C-1, C-2 (set C) and P3, C-2 (H-pool) next to the pre-specified superfamily-weighted ones, from stored per-pair rows | `results/r44_pairw.json` |
| `r44_p1win.py` | Thermophile statistic U2 (Q ≥ 0.7) in absolute 320 K windows of 50–100, 300–400 and 400–500 ns on the same pairs and shared bonds as P1′ (spread → stats) | `derived/r44/p1win/`, `results/r44_p1win.json` |
| `r44_corr2.py` | Agreement of the correspondences against TM-align on medoids, an alignment independent of π̄, with tolerances of 0–2 bonds by identity (r44_agree); TM-align on the α-carbon traces of the projected mean points and C-1, C-2 on that correspondence (r44_tmmean) | `derived/r44/corr2/pairs_c2_s*.jsonl`, `results/r44_corr2.json` |
| `r44_b6size_cell26.py` | Partial Spearman of Δ²_D/√n̄ with TM-score given log n̄ on the continuous axis (r44_b6size). Run in a notebook cell; the file keeps the cell verbatim below a short header that loads the modules | `results/r44_b6size.json` |
| `r44_seqid_cell30.py` | Implementation check after the r44_corr result (not specified in advance, not used in the manuscript): share of π̄ bond pairs contained in the sequence alignment by π̄ identity. Notebook cell kept verbatim | standard output |
| `r44_repro.py` | Synthetic check that framedchain 0.4.1 reproduces the as-run code (framespace 0.2.0, `framespace_wip/stage0`, `red_pairs` and `estats` of `r38_s0.py`) exactly: distances, alignments, certificates, mean points, ensemble W₂, reduced computation, energy statistics. Run as `python code/r44_repro.py <dir>` where `<dir>` holds checkouts of this repository (`framedchain-paper/`) and of the package (`framedchain/`) side by side | standard output |
| `r44_w2p1.py` | Part w2: why W₂ and Δ²_D rank the 122 continuous-axis pairs differently — Spearman correlations with the replicate spread w̄ and the cross-ensemble mean distance C, and the partial correlation of W₂ with Δ²_D given w̄ (r44_w2gap). Its part p1comp stopped before producing values (a field it read was never filled) and was replaced by `r44_p1comp2.py` | `results/r44_w2gap.json` |
| `r44_p1comp2.py` | Composition of thermophile–mesophile and control pairs in the P1′ replication (log shared bonds, shared fraction, normalized edit distance between mean points, length), and U2 after removing a linear dependence on them fitted on the control pairs (r44_p1comp, as corrected before it was computed) | `results/r44_p1comp.json` |
| `r44_redsame.py` | Whether the reduced computation produces the compression of Δ²_D relative to aligned coordinates on same-sequence cells: in 16 E7 domains, 20 frame pairs each between 320 and 379 K and within each temperature, compared with branch-and-bound certified values; share of the compression attributable to the reduced computation (r44_redsame) | `results/r44_redsame.json` |
| `r45_f1.py` | Whether the sequence-identity effect of set C survives when the simulated fluctuation profiles are restricted to bond directions: from the same GPA superposition, whole-frame spreads (as pre-specified), the spread of the unit vector f₁ (tilt) and the rotation about the bond [MSD(f₂) + MSD(f₃) − MSD(f₁)]/2 (twist); C-1 and C-2 for the three versions with the same superfamily bootstrap weights and paired differences, after an implementation check that the whole-frame version reproduces the recorded values; agreement of the network-model predictions with each version (corr K k → sanity → stats) | `derived/r45/f1/`, `results/r45_f1.json` |
| `r45_mpstart.py` | Independence of the mean point from the starting frame with the as-run `fs.mean_point`: six set-C ensembles rebuilt with `r38_c.dom_build`, nine starting frames each, distances after one optimal rotation computed directly, and the start-0 mean point compared with the stored one | `results/r45_mpstart.json` |
| `r45_folded_frac.py` | Descriptive count (not specified in advance) for the 67 P2 domains: share of MD 320 K frames over whole trajectories with Q ≥ 0.75 (the `md_load` definition) and number of folded aSAMt samples. Run from a notebook as `/content/folded_frac.py`; the file keeps it verbatim below a header added for this repository | `results/r45_folded_frac.json` |

## Round46–49 — checks added after external review (2026-10-29 to 10-31; result 4.5 and Methods)

Each analysis was written into a new time-stamped record before it was run (`registration/prereg_2026-10-30_addendum*.json`
and `prereg_2026-10-31_addendum12*.json`); corrections were written as separate records before the affected step was computed.

| Script | What it does | Main outputs |
|---|---|---|
| `r46_align.py` | Whether alignment accuracy shared by sequence identity and r_dyn accounts for their association in set C. Re-runs TM-align on the medoid α-carbons (tm) and checks it against `r44_corr2` integer by integer; then (stats) A = identity over all residues aligned by TM-align, B = the precision of π̄ against TM-align (exact, within two bonds) added to the C-1 and C-2 controls, with descriptive versions, from one superfamily bootstrap stream of 10,000 (seed 20261121). Its second check failed because tied identity values were ordered differently by `argsort` on different runtimes | `derived/r46/align_pairs.jsonl`, `derived/r46/tm_check.json`, `results/r46_align.json` |
| `r46_align_fix1.py` | The statistics of `r46_align.py` with average ranks (`scipy.stats.rankdata`), as fixed by the correction record; it loads the functions of `r46_align.py` and replaces only the rank function. Checks 2a (argsort version against `results/r38_c.json`) and 2b (average-rank baseline against the values written in the correction) | `results/r46_align_fix1.json` |
| `r47_ties.py` | Effect of tie handling on rank-based estimates whose intervals end near zero: the tie structure of the variables (ties), and the first homolog study with its superfamily-weight interval (p3), the linear shorter-domain controls of `r38_c_nl` (nl) and the network-model versions of `r43_sens` s3 (anm), each with the original ranks (argsort) and with average ranks; combine compares both with the stored values. The parts run in parallel from `run_records/queues/q_r47.py` | `derived/r47/<part>_<mode>.json`, `results/r47_ties.json` |
| `r47_combine_fix1.py` | Comparison step of `r47_ties.py` reading only estimate fields (the original step had read text fields as estimates). It stopped with a syntax error because it cut `r47_ties.py` at the first `if __name__`, which occurs inside a string | none |
| `r47_combine_fix2.py` | The same comparison step, cutting `r47_ties.py` at the last `if __name__` | `results/r47_ties.json` |
| `r48_identity_range.py` | Whether the association in set C extends below 15% identity: computes the composition-expected identity Σ_a p_a q_a from the domain files (exp), then C-1 and C-2 below 15% (decision set) and 10%, with identity from TM-align, with the precision covariate, over all pairs of the superfamilies that have a pair below 15%, values expected under range restriction (1,000 simulated data sets), residual r_dyn by identity decile and observed against expected identity (stats); average ranks, one bootstrap stream of 10,000 (seed 20261131) | `derived/r48/exp_ident.jsonl`, `results/r48_identity_range.json` |
| `r49_chainbreak.py` | Chain-break check (Cα–Cα > 4.2 Å) of the sets not checked on 2026-10-28, reading the atoms each analysis read: the mdCATH files (first frame of each 320 K replicate; for the thermophile replication, the stored second-half frames where present), the shared residues of the ATLAS starting structures, and the static chains of the first two Results sections. Its run stopped when the Colab VM was replaced | none |
| `r49_chainbreak_fix1.py` | The same computation, writing each domain or chain to `derived/r49/items.jsonl` so that a rerun skips finished items, with 16 workers | `results/r49_chainbreak.json` |

## Figures and tables

These scripts read stored results only (see `docs/RESULTS_MAP.md`, Figures and tables); on Colab, the figures (in Arial when installed, with msttcorefonts, otherwise Liberation Sans; text 8 pt or larger) and the tables were last drawn as follows: all eight figures together on 2026-11-02 at 06:07 UTC, the S1 and S3 Tables on 2026-11-02 at 06:07 UTC and the S2 Table on 2026-11-01 at 18:29 UTC (written again with the same content on 2026-11-02 at 06:07 UTC; times of the files on Colab; earlier runs in `docs/PROVENANCE.md`).
They take no arguments; the project path `R` and the output folder are set at the top of each file.

| Script | What it does | Main outputs |
|---|---|---|
| `make_figs.py` | Figures 1–7 and S1 Fig: RGB TIFF with LZW compression at 600 dpi and preview PNG files | `manuscript_figs/Fig1.tif` … `Fig7.tif`, `S1_Fig.tif` |
| `make_tables.py` | S2 Table (κ and Δ²_D/√n̄ of all and size-matched cells) and S3 Table (alignment-accuracy and identity-range checks) | `manuscript_figs/S2_Table.xlsx`, `S3_Table.xlsx` (+ .csv) |
| `make_s1.py` | S1 Table: entries of the first record with their times (sheet A), the added records (B), the pre-specified predictions with outcomes (C) and the classification of every analysis in the manuscript (D), with notes | `manuscript_figs/S1_Table.xlsx` |

## Node.js checks

| Script | What it does |
|---|---|
| `negtype/negtype.js` | Tests whether D is of negative type on small random structures (positive eigenvalue of J D J means it is not) |
| `negtype/so3fourier.js` | For length-2 structures D depends only on the angle θ between relative rotations; computes the SO(3) Fourier coefficients of f(θ) = h·min(√6, 4 sin(θ/4)) (a positive coefficient rules out negative type) |
| `negtype/negtype_eta_len2.js` | η₊ (positive spectral mass of J D J) for 60, 100, 200 and 400 random two-bond chains, on which D is provably not of negative type, using the functions of `negtype.js`; random rotations from the fixed-seed generator of `negtype.js` (seed 12345), so every run prints the same values (S3 Text 2). Run as `node negtype_eta_len2.js .` from `code/negtype` |
| `negtype/len2check.js` | Checks the closed form of D for length-2 structures against the general computation (maximum difference about 1e-15), then prints the top eigenvalues of J D J for 60, 120 and 200 random length-2 structures (positive values appear from 120 on) |

Run them from `code/negtype` (`node so3fourier.js`, `node len2check.js`, `node negtype.js`); the printed labels
are in Korean. `so3fourier.js` prints the coefficients for D (with the length cap), without the cap, and for
the chordal distance of rotation matrices as a control; the positive coefficient f̂_5 = +1.15 × 10⁻² is the
value cited in result 1.1.
