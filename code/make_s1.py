"""S1 Table: time-stamped records, pre-specified predictions with outcomes, and the classification of every analysis in the paper."""
import json, glob, os, re, pandas as pd
R = '/content/drive/MyDrive/atlas_framed_curves'; OUT = R + '/manuscript_figs'
USED = {'round17_chaining', 'round29_review_math', 'round37', 'round42_withdrawn_2026-10-21', 'round43_2026-10-22', 'round41_2026-10-21'}
FLAG = {'round36_p1_checklist_2026-10-20': 'Flagged: stamped 12:57:26 inside the entry, but the entry reached Drive at 12:57:39, 3 s after the decision file r38_p1.json (12:57:36). The placeholder "running" for the condition-2 simulation (its result file created 12:57:26) was replaced by its summary at 12:58:12 (stamp inside 12:58:00, noted as before reading r38_p1.json), 24 s after the decision file by that stamp and 37 s by Drive time; because the text placeholder became a list, the change is listed under Items removed.',
        'round36_roundtrip_2026-10-20': 'Flagged: the entry reached Drive (17:34:46) after its output r39_roundtrip.json (17:26:04); the stamp inside the entry is 17:23:15 and the prediction is an identity (maximum round-trip error < 1e-9 Å).',
        'round41_2026-10-21': 'Two checksums of a data-preparation script edited in place (07:44 and 07:52 UTC on 22 October), before its outputs were produced (08:42).',
        'round43_2026-10-22': 'Checksum and stamp of the post hoc sensitivity script edited in place 70 s after its output file was created (08:59:18 vs 08:58:08).',
        'round34': 'A partial outcome replaced by its completed version (16:03:55).',
        'round36_partial_freeze_7b7c': 'Written after the 7b split outputs (1,640 files, 19 Oct 11:13–14:53) were computed and before they were read (state recorded in the entry); predictions 1–3 of 7b were in the record from 18 Oct 18:58:51.',
        'round36_p3_spec_2026-10-19': 'A smoke test of r38_p3.py computed values for one pair to check that they were finite; its log was saved at 18:37:52, after the entry reached Drive (snapshot copy 18:37:38, retained version 18:37:46). The script and its decision rules were fixed before (Drive copy last modified 18:37:02).',
        'round36_trackB_code_2026-10-19': 'First entered at 19:42:45 (retained Drive version; stamp inside 19:42:33) under the key round36_trackB_code_2026-10-20, with the same script checksums (r38_b6.py among them) and rules; re-stamped under the present key at 19:43:30 with only the date in the key changed to UTC (the entry records this in note and first_stamp_utc). First seen on Drive refers to the present key. Earlier files: r38_verify_b.json and r38_verify_b_run2.json (18:54–19:38) are implementation checks of r38_p3.py and r38_b6.py on synthetic inputs and on single real domains and pairs, without decision statistics; r38_s0_bias2_smoke.json (19:23:11) is a code test of the re-measurement that timed two settings and reproduced the original bias from stored errors, without computing the bias of any new setting (recorded as NaN); the choice rule was already in the header of r38_s0_bias2.py.',
        'round36_p1_code_2026-10-20': 'The pair-counting step (04:56:29) and the start of the mean-point step (04:58:54) of r38_p1.py ran before the entry reached Drive (04:58:59; stamp inside 04:58:52); the decision statistics ran at 12:57. r38_verify_p1.json (04:56:52) holds implementation checks of r38_p1.py: nine on synthetic inputs and three on real data (frame counts for one domain, alignment for one pair, and the full computation for one pair, for which only the type of θ*, not its value, was stored).',
        'round36_amend_C_2026-10-20': 'Data download for the confirmatory set began at 06:17, before the entry (state at stamp: 42 of 883 files); 7g (E7 on 100 domains) was specified in the entry before its run (06:44), and its script checksums were added at 06:50:51, after the run had started and before its result (15:33). r38_c_probe.json (06:16:05) is a data-size probe (file structure, number of domains to download, estimated time) without statistics.',
        'round36_p1_diag_2026-10-20': 'Stamped 13:15:29 inside; the diagnostic ran from 13:15:32, before the entry reached Drive at 13:15:38, and after the P1 result (12:57:36) it diagnoses.',
        'round39_2026-10-20': 'Stamped 17:18:31 inside; reached Drive at 17:34:46 together with three other entries. A preparation step for E5 (aligned coordinates, no statistics) ran from 17:19:22; the E0′–E5 statistics ran from 17:40:28.',
        'round36_sec6_2026-10-20': 'Reached Drive at 17:34:46; the 128-start alignment step of the static comparison began at 17:33:43 (machine clock) and branch-and-bound at 17:55:07; the predictions concern the certification.'}
INFILE = {'round34': '2026-10-18 (date only, KST; written locally before running, appended unchanged)', 'round35': '2026-10-19 (date only, KST; before running)',
          'round37': '2026-10-19 (date only, KST; before running)', 'round36': '2026-10-19 (date only, KST; third version)'}
def entries():
    t = pd.read_csv(f'{R}/results/_registration_evidence_2026-10-26/entries_timeline.tsv', sep='\t').fillna('')
    t = t[~t.entry.isin(['written', 'note'])].copy()
    t['Bears on this paper'] = ['yes' if (e in USED or e.startswith('round36') or e.startswith('round39')) else 'no (earlier exploratory stage)' for e in t.entry]
    t['Note'] = [FLAG.get(e, '') for e in t.entry]
    t['in_file_time'] = [INFILE.get(e, v if v else '— (record begun 2026-10-17 12:47:55)') for e, v in zip(t.entry, t.in_file_time)]
    t = t.rename(columns={'entry': 'Entry', 'first_seen_utc': 'First seen on Drive (UTC)', 'n_versions': 'Versions', 'change_times_utc': 'Later changes (UTC)',
                          'in_place_modified_paths': 'Values edited in place', 'removed_paths': 'Items removed', 'result_files_created_before_first_seen': 'Result files created earlier (Drive time)', 'in_file_time': 'Time written inside the entry (UTC; date-only values in KST)'})
    return t[['Entry', 'Time written inside the entry (UTC; date-only values in KST)', 'First seen on Drive (UTC)', 'Versions', 'Later changes (UTC)', 'Values edited in place', 'Items removed', 'Result files created earlier (Drive time)', 'Bears on this paper', 'Note']]
ADDED = [  # file, written (UTC), analyses, what it fixed
 ('prereg_2026-10-28_addendum.json', None, 'r44_inv, r44_corr', 'Replicate- and frame-number audits (±10% criterion) and the decomposition of κ² into mean-shift and spread terms; alternative correspondences (TM-align, sequence alignments).'),
 ('prereg_2026-10-28_addendum2.json', None, 'r44_corr_v2, r44_gp, r44_pairw, r44_p1win', 'Correspondence run split into eight parallel chunks to include the glycine/proline version (computation, statistics and decisions unchanged; the first run stopped unseen); glycine and proline; pair-weighted estimates; absolute trajectory windows.'),
 ('prereg_2026-10-28_addendum3.json', None, 'r44_b6size', 'Size control in the comparison across sequences.'),
 ('prereg_2026-10-28_addendum4.json', None, 'r44_agree, r44_tmmean', 'Agreement of each correspondence with TM-align on medoids; TM-align on projected mean points.'),
 ('prereg_2026-10-28_addendum5.json', None, 'r44_w2gap, r44_p1comp', 'Source of the gap between W2 and Δ²_D; composition of thermophile and control pairs.'),
 ('prereg_2026-10-28_addendum5_fix1.json', None, 'r44_p1comp (correction)', 'Covariates TM-score and identity were not recorded for these pairs; replaced by D̂ between mean points, before the analysis was computed.'),
 ('prereg_2026-10-28_addendum6.json', None, 'r44_redsame', 'Contribution of the reduced computation to the compression of d on same-sequence cells.'),
 ('prereg_2026-10-29_addendum7.json', None, 'r45_f1', 'Fluctuation profiles restricted to bond directions; rotation about the bond.'),
 ('prereg_2026-10-29_addendum8.json', None, 'r45_mpstart', 'Dependence of mean points on the starting frame.'),
 ('prereg_2026-10-30_addendum9.json', None, 'r46_align', 'Whether alignment accuracy shared by identity and profile similarity accounts for their association.'),
 ('prereg_2026-10-30_addendum9_fix1.json', None, 'r46_align (correction)', 'After the check of recorded estimates failed (tie order differs between processors) and before comparisons: redefined with average ranks.'),
 ('prereg_2026-10-30_addendum10.json', None, 'r47_ties', 'Effect of tie handling on rank-based estimates with intervals ending within 0.011 of zero.'),
 ('prereg_2026-10-30_addendum10_fix1.json', None, 'r47_ties (correction 1)', 'Comparison step had read text fields as estimates; computations unchanged.'),
 ('prereg_2026-10-30_addendum10_fix2.json', None, 'r47_ties (correction 2)', 'The first correction failed to load the script named in the record; computations unchanged.'),
 ('prereg_2026-10-30_addendum11.json', None, 'r48_identity_range', 'Whether the association extends below 15% identity, with sentences for each outcome.'),
 ('prereg_2026-10-31_addendum12.json', None, 'r49_chainbreak', 'Chain-break check (Cα–Cα > 4.2 Å) of the sets not checked on 28 October.'),
 ('prereg_2026-10-31_addendum12_fix1.json', None, 'r49_chainbreak (correction)', 'After the first run stopped without output when the Colab virtual machine was replaced: the same computation, resumable item by item, with 16 workers.'),
]
def added():
    rows = []
    for f, _, an, what in ADDED:
        A = json.load(open(f'{R}/results/{f}')); w = A.get('written') or next((v.get('written') for v in A.values() if isinstance(v, dict) and v.get('written')), '')
        rows.append({'File (date in KST)': f, 'Written (UTC, inside the file)': w, 'Analyses': an, 'What it fixed': what})
    return pd.DataFrame(rows)
L_DC, L_SC, L_D, L_BD = 'decision criterion', 'sanity check', 'descriptive', 'by definition'
PRED = [  # record, analysis, prediction, label, outcome, where reported
 ('round17_chaining (18 Oct 00:20 UTC)', 'Chaining test, 150 two-domain chains', '1 − TM-score normalized by the shorter chain violates the triangle inequality in > 80% of part–whole triples; D_TM with fixed d₀ and D_TM/√L̄ near 0.', L_DC, 'Held: 100% vs 0% and 0%.', 'R2, Fig 2B'),
 ('round17_chaining', 'Chaining test', 'Average linkage merges different folds under 1 − TM-score (shorter chain) but mostly not under D_TM/√L̄; single linkage merges under both.', L_DC, 'Held: average linkage 100% vs 15% [10, 21]; single linkage merged under both (100% and 78%), as predicted.', 'R2, Fig 2B'),
 ('round29_review_math (18 Oct 08:40 UTC)', 'Chaining test, fixed d₀', 'D_TM/√L̄ with a d₀ common to each triple: no violations, average-linkage chaining as for the d₀(L̄) version (about 15%).', L_DC, 'Held: 0% and 15.3% [10.0, 21.3].', 'R2, Fig 2B'),
 ('round29_review_math', 'Chaining test, linear form', 'Linear 1 − S/L̄ violates the triangle inequality and chains more than the square-root form.', L_DC, 'Partly: violations 66.7% (held); chaining 9.3% (not held).', 'R2, Fig 2B'),
 ('round36 §0 (full text in record 19 Oct 14:56 UTC)', 'First session 0.0', 'Branch-and-bound and heuristic agree on test pairs.', L_SC, 'Passed (BnB 10/10).', 'S2 Text'),
 ('round36 §0', 'First session 0.1 (truncated chains)', 'Unshared part of full-versus-truncated minus full-versus-full comparisons equal to 0.75 Σ s² over the removed bonds within 2% in every domain, and median shared-part Δ²/√n within the replicate band.', L_DC, 'Not held (relative errors −6% to −107%; terminal flexibility).', 'S2 Text'),
 ('round36 §0', 'First session 0.2 (upper bias of the reduced computation)', 'Bias below 5% of Δ²_D in every TM-score bin.', L_DC, 'One bin at 1.02× the threshold; band widened under the pre-specified rule, with settings fixed after the result and before re-measurement (0.70× after).', 'Methods, S2 Text'),
 ('round36 §0', 'First session 0.3 and 0.4', 'Incomplete U-statistic agrees with the complete statistic; no flags.', L_SC, 'Passed; flags 0/28.', 'Methods, S2 Text'),
 ('round36 §0, §3 (full text in record 19 Oct 14:56 UTC)', 'E7 (5 and 30 domains)', 'On same-sequence cells, Spearman ≥ 0.9 between the edit-distance and aligned-coordinate versions of Δ², w̄ and the ratio, and median ratio of Δ²_D to aligned-coordinate Δ² within [0.8, 1.2].', L_DC, 'Not held: 0.65 and 0.72 (ranks preserved, ρ = 0.98); pre-specified rule: magnitudes with d.', 'R5'),
 ('round36 §1 (rules in r38_b6.py, checksum in record 19 Oct 19:42 UTC)', 'Continuous axis, prediction 1', 'Spearman(Δ²_D/√n̄, TM-score) < −0.5.', L_SC, 'Held: −0.83 [−0.89, −0.75] (bootstrap over domains; the decision used the superfamily-cluster interval [−0.87, −0.78]).', 'R4'),
 ('round36 §1', 'Continuous axis, prediction 2', 'Δ²_D/√n̄ decreases across TM-score bins without saturating.', L_D, 'Held: 0.83 to 0.23.', 'R4'),
 ('round36 §1', 'Continuous axis, prediction 3', 'Far-limit prediction relative residual < 15% for TM-score < 0.5 (larger for TM-score ≥ 0.7).', L_DC, 'Held: 3.7% (19%).', 'R4'),
 ('round36 §1', 'Continuous axis, prediction 4', 'Unshared part per bond nearly independent of TM-score, shared part per bond larger at lower TM-score (descriptive, because of the shrinkage bias of W₂²).', L_D, 'Unshared part not independent of TM-score (ρ −0.87; not held); shared part larger at lower TM-score (ρ −0.67; held).', 'not discussed'),
 ('round36 §1', 'Continuous axis, prediction 5', 'η₊ ≤ 0.10 in at least 80% of same-superfamily pairs.', L_DC, 'Held: 100% (also in all 234 pairs).', 'R4'),
 ('round36 §1', 'Continuous axis, prediction 6', 'Spearman(Δ²_D, W₂) ≥ 0.8 for TM-score ≥ 0.5.', L_DC, 'Not held: 0.74; W₂ shown alongside Δ²_D (pre-specified rule).', 'R4'),
 ('round36_r10a (20 Oct 16:55 UTC)', 'Natural insertion–deletion cell (R10a′), prediction 3', 'Shared part per bond of indel pairs lies on the homolog regression line (equivalence margin ±0.0146).', L_DC, 'Not held: residual −0.089 [−0.105, −0.043] ("different"); not separable from sample differences after controls (added analysis).', 'not discussed (cell shown in Fig 4B)'),
 ('round36 §2; criteria in round36_r10a (20 Oct 16:55 UTC)', 'R10a′, prediction 2', 'Unshared share minus self-comparison equals 0.75 × inserted residues (relative error < 10%).', L_DC, 'Not held (median relative error 24-fold).', 'not discussed'),
 ('round36_sec6 (20 Oct 17:34 UTC)', 'Certification, prediction 1', 'At least 95% of 500 pairs certified within 60 s.', L_DC, 'Held: 500/500.', 'R1'),
 ('round36 §6', 'Certification, prediction 2', '128-start heuristic within 0.002 of the certified value in ≥ 95% of pairs, as (S_ub − S_heur)/S_ub ≤ 0.002.', L_DC, 'Specification error: could not hold by construction (S_ub = S_lb + tolerance); 0/500 as specified, 499/500 within 0.002 of the best value found (post hoc).', 'R1 (values), S2 Text'),
 ('round36 §7 (rules in r38_p1.py, checksum in record 20 Oct 04:58 UTC)', 'Thermophile P1, condition 1', 'Superfamily-median θ* with one-sided 97.5% lower bound > 320 K.', L_DC, 'Not held: 333 K, lower bound < 320 K (40% left-censored).', 'R7'),
 ('round36 §7', 'Thermophile P1, condition 2', 'Paired difference from mesophile–mesophile control, lower bound > 0.', L_DC, 'Held: +21.2 K (lower bound +4.3 K).', 'R7'),
 ('round36 7b (full text in record 18 Oct 18:58 UTC)', 'Standard metrics on public AlphaFlow-MD samples, prediction 1', 'Energy distance calls every controlled perturbation far (≥ 95%).', L_DC, 'Not held: inflation 0.93, truncation 0.62 and 0.70 (others 0.99–1.00).', 'not discussed'),
 ('round36 7b', '7b, predictions 2 and 3', 'RMWD calls shrinkage and truncation far (≥ 70%); contact Jaccard calls mixing far (≥ 90%).', L_DC, 'Held: 0.81–0.99; 1.00.', 'not discussed'),
 ('round36 7c (full text in record 18 Oct 18:58 UTC)', 'Model error P2 (primary)', 'Heating zone in which the superfamily median of the aSAMt mean-shift error falls.', L_DC, 'Censored zone "≥ +59 K".', 'R6, Fig 5A'),
 ('round36 7c', 'Model error, secondary 1', 'Error above replicate noise in at least 70% of domains.', L_DC, 'Held narrowly: 72.7% (interval [61.5, 83.3] from an added check).', 'R6'),
 ('round36 7c; edit-distance version in round36_amend_P2_D2 (20 Oct 11:25 UTC)', 'Model error, secondaries 2 and 3', 'AlphaFlow mean shift exceeds the force-field threshold; position of the aSAMt error in the force-field bracket.', L_DC, 'Held (0.0843 > 0.0199); "outside the bracket" (p = −0.091).', 'not discussed'),
 ('round36 7c', 'Direction of the error (with round39 E5)', 'Cosine and profile correlation judged against ceiling and replicate null.', L_DC, 'Mean shift: different direction (0.25); profile: same direction at the highest temperature, different at 379 K.', 'R6, Fig 5B'),
 ('round36 7h (20 Oct 13:36 UTC)', 'Thermophile replication P1′-1', 'U2 (Q ≥ 0.70) one-sided 97.5% upper bound < 0.', L_DC, 'Held: −0.213 (upper bound −0.083).', 'R7, Fig 6'),
 ('round36 7h', 'P1′-2 (decomposition)', 'U2(Q ≥ 0.70) − U2(Q ≥ 0.80) upper bound < 0.', L_DC, 'Not held: −0.016 [−0.144, 0.019].', 'R7'),
 ('round36 7h', 'P1′-3 (stability)', 'Paired difference in folded-fraction logit at 379 K, lower bound > 0.', L_DC, 'Held: +1.58 (lower bound 1.05).', 'R7, Fig 6'),
 ('round36 1b (rules in r38_p3.py, checksum in record 19 Oct 18:37 UTC)', 'Sequence identity P3 (first homolog study)', 'Partial Spearman r_dyn ~ identity | shared-core structure. Effect: lower bound > 0, upper bound ≥ 0.15 and lower bound > 0 also with the intersection of near-optimal alignments on the same pairs; smaller than the minimum effect: lower bound > 0, upper bound < 0.15; equivalent: interval within ±0.15; otherwise inconclusive (equivalence within ±0.10 reported separately).', L_DC, 'Inconclusive: 0.111 [−0.048, 0.266].', 'R8, Fig 7'),
 ('round36 1c (20 Oct 06:23 UTC)', 'Confirmatory C-1 and C-2', 'Lower bounds > 0 with the structural control and with the network-model control added.', L_DC, 'Held: 0.224 [0.119, 0.333]; 0.181 [0.079, 0.283].', 'R8, Fig 7'),
 ('round36 1c', 'Auxiliary analyses', 'Superfamily size; log length as covariate.', L_D, '0.291 (≤ 2 pairs) vs 0.167; 0.218 with log length.', 'R8'),
 ('round36_c_nl (20 Oct 15:26 UTC)', 'Nonlinear controls (added)', 'After cross-fitted gradient boosting on eight covariates: lower bound > 0 → the effect remains; lower bound ≤ 0 and upper bound < 0.15 → explained by nonlinear functions of the structural variables; otherwise inconclusive.', L_DC, 'Held: 0.141 [0.045, 0.228]; shorter domains 0.032 [−0.049, 0.109].', 'R8, Fig 7'),
 ('round36_roundtrip', 'Coordinate round trip', 'Maximum Cα round-trip error < 1e-9 Å.', L_BD, 'Held: 8.0e-15 Å.', 'R1'),
 ('round36_r6', 'Three same-protein events (round36_r6)', 'Rank correlation of shared part with wild-type spread ≤ −0.5.', L_D, 'ρ = −0.5 (meets ρ ≤ −0.5; three events).', 'not discussed'),
 ('round37 (18 Oct 17:08 UTC)', 'round37 R3″ (protocol only)', 'Median Δ² of protocol-only differences at most half of the force-field-plus-protocol value (0.076), that is ≤ 0.038 (reported at a descriptive level).', L_D, 'Held (within the replicate-only band).', 'not discussed'),
 ('round37', 'round37 R3′ (force field and protocol, same protein)', 'Temperature equivalent of the mean shift at least +59 K after subtracting 20 K, reported regardless of direction.', L_D, 'Whole trajectories +90 K (held); second halves +58 K (just below); not converted to temperature in the paper because the direction differs from heating (E5).', 'not discussed (R4 reports κ)'),
 ('round37', 'round37 R5′ (substitutions)', 'Temperature equivalents of single substitutions between +28 and +59 K (reported; force fields differ).', L_D, 'Held in 1 of 4 systems.', 'not discussed'),
 ('round37', 'round37 R7 (upper anchor)', 'Folded versus unfolded Δ²/√n at least 5 times the +130 K cell (0.066); its placement on the scale to be decided later.', L_D, 'Split (Trp-cage only in one set; both in the mdCATH set).', 'R3 (no stable upper anchor)'),
 ('round37', 'Unit invariance audit', 'Second-half / whole and survivor / all ratios within [0.8, 1.25]; otherwise the second half (or the restricted set) becomes primary.', L_DC, 'Window not met (1.29, 1.33 and 1.01 at +59, +93 and +130 K; 2.10 at +28 K; looser folded-frame criterion; second halves used as primary); survival mostly met (1.05 and 0.95 at +59 and +93 K; 1.59 at +28 K, outside the range).', 'Table 2'),
 ('round39 (20 Oct 17:34 UTC)', 'E1 size exponent', 'α ordered 59 < 93 < 130 K, α(+130) − α(+59) lower bound > 0, intervals exclude 0.5 (all three needed).', L_DC, 'Not held (order held; difference 0.15 [−0.57, 0.82]; 0.5 excluded).', 'R3, Fig 3B'),
 ('round39', 'E3 resolution, predictions 1 and 2', 'At 20 domains, +59 K and each larger step detected with probability ≥ 0.8; +28 K needs at least 100 domains for 0.8.', L_DC, 'Prediction 1 not held (0.38, 0.73 and 0.92 at +59, +93 and +130 K); prediction 2 held (200 domains).', 'R3, Fig 3C'),
 ('round39', 'E4 gate', 'Gate: superfamily-median temperature equivalents at 379 and 413 K within 10 K; mean absolute error ordered shape × amplitude ≤ linear interpolation < population curve; 450 K version within 10 K.', L_DC, 'Gate not met in second halves, met over whole trajectories; 450 K not met; error ordering not held; equivalents reliable up to +93 K.', 'R3'),
 ('round39', 'E5 force-field cell (sanity check)', 'Spread change of the force-field cell opposite to heating.', L_SC, 'Held: −0.54 [−0.63, −0.36].', 'not discussed'),
 ('round39_amend1 (21 Oct 08:27 UTC)', 'E6 upper anchor', 'Folded-versus-unfolded cells give a stable upper anchor.', L_DC, 'No anchor (τ = 0.943 held, position stability not held).', 'R3'),
 ('round39', 'E0′ regimes, E2 window', 'Regime of each cell; window ratios.', L_D, 'Reported.', 'R3, R4, Fig 4B, Table 2'),
 ('round41 (21 Oct 08:31 UTC)', 'g1 static axis S2 and S3', 'S2: partial Spearman of the unshared part per √n̄ with the CATH level separating the pair (H < T < A < C < X), given TM-score and length difference, lower bound > 0. S3: adding the unshared part to TM-score and length difference improves the classification of same-superfamily pairs, ΔAUC lower bound > 0.', L_DC, 'Not held: S2 −0.024 [−0.084, 0.067]; ΔAUC +0.0004 [−0.0002, 0.0010] (AUC with TM-score alone 0.977).', 'not discussed'),
 ('round41', 'g1 unshared share and r_dyn', 'The sign of the association between the unshared part of the mean-point alignment (per √n̄) and r_dyn, controlling for medoid TM-score and length difference, found in the first homolog set, is confirmed in the confirmatory set (95% interval excluding 0 in the same direction).', L_DC, 'Confirmed: −0.146 [−0.241, −0.049]; with the shared-part control as well, −0.078 [−0.167, 0.009].', 'not discussed'),
 ('round41', 'g2 compression of d (added after the E7 result)', 'Reading rule: the step carrying more than half of the log compression names the mechanism (rotation per pair, saturation beyond 120°, or shifts of the correspondence).', L_D, 'Shifts of the correspondence 0.72 [0.68, 0.75]; saturation 0.29; rotation −0.015.', 'R5'),
 ('round41', 'g6 synthetic sequences', 'No violations for 1-Lipschitz gap costs (i)–(iii) (0 of 3,000 inequalities); violations expected for the controls (iv) 2|a| and (v) DTW, without refuting the lemma if absent.', L_DC, 'Held for (i)–(iii); the controls also had no violations under random edits.', 'S1 Fig (adversarial version)'),
 ('round41 g6_adv (21 Oct 09:25 UTC; added after g6)', 'g6 adversarial triples', '(i)–(iii) no violations; (iv) violations on generator A; (v) violations on generator B.', L_DC, 'Held: 0; 21.1% and 22.7%; 33.0% (65 undecided equalities in (iii) A).', 'R1, S1 Fig'),
 ('round42 (21 Oct 09:27 UTC)', 'Hinge-type test', 'Withdrawn before running (premise contradicted by a decomposition of explanatory variables).', '—', 'Withdrawn.', 'not discussed'),
 ('round43 (22 Oct 08:54 UTC)', 'Post hoc sensitivities s1, s3, s4', 's1: if a cell changes regime under length matching or Kendall τ < 0.9, show both versions; s3: if C-2 and its intersection gate hold for all four network models, C-2 is not tied to the parameters; s4: interval of secondary 1, verdict unchanged.', L_D, 's1: no regime changed in the 13 cells with size-matched units, τ = 0.91 over 12 (matched version in S2 Table); s3: effect in all four (C-2 0.18–0.19); s4: [61.5, 83.3], includes 0.70.', 'R4, R6, R8, S2 Table'),
 ('addendum (28 Oct 04:47 UTC)', 'r44_inv', 'Reduced/original ratio of superfamily medians in [0.9, 1.1] with 95% interval including 1, at +59 and +93 K, for 3 vs 5 replicates and 50 vs 100 frames.', L_DC, 'Replicate axis not met (+93 K 0.87 [0.45, 1.22]); frame axis met.', 'R3, Table 2'),
 ('addendum', 'r44_corr', 'P1: TM-align contains π̄ pairs more than sequence alignments; P2: TM-align C-1, C-2 lower bounds > 0; P3: sequence alignments lower bounds > 0 with C-1 below π̄.', L_DC, 'P1 held; P2 held; P3 not held (C-1 below π̄ held).', 'R8, Table 3'),
 ('addendum2 (28 Oct 05:56 UTC)', 'r44_gp', 'Effect remains without bonds next to Gly/Pro and with Gly/Pro fractions as covariates (lower bounds > 0, more than half of C-1).', L_DC, 'Held.', 'R8, Table 3'),
 ('addendum2', 'r44_pairw', 'P1: pair-weighted C-1 below superfamily-weighted; P2: its lower bound > 0; P3: pair-weighted interval of the first-study test P3 includes 0.', L_DC, 'P1 held; P2 held; P3 not held (0.077 [0.043, 0.100]).', 'R8 (P1, P2); P3 not discussed'),
 ('addendum2', 'r44_p1win', 'PA1: U2 < 0 in both late windows; PA2: 300–400 ns upper bound < 0.', L_DC, 'Held; held.', 'R7, Fig 6'),
 ('addendum3 (28 Oct 06:04 UTC)', 'r44_b6size', 'Partial Spearman controlling for log n̄ ≤ −0.6 with upper bound < 0.', L_DC, 'Held: −0.84 [−0.89, −0.76].', 'R4'),
 ('addendum4 (28 Oct 06:46 UTC)', 'r44_agree', 'PA: below 15% identity, sequence alignments contain < 0.2 of TM-align bond pairs (median, ±2); PB: π̄ contains more than sequence alignments.', L_DC, 'Held (0.00); held (0.66).', 'R8, Table 3'),
 ('addendum4', 'r44_tmmean', 'Reading rule for the share of the C-1 difference recovered by aligning mean points.', L_D, 'Share 0.58: "both" (alignment method and structures aligned).', 'R8'),
 ('addendum5 (28 Oct 08:53 UTC)', 'r44_w2gap', 'PW1: partial Spearman(W₂, Δ²_D | w̄) ≥ 0.8; PW2: ρ(W₂, w̄) > ρ(Δ²_D, w̄).', L_DC, 'Held (0.95); held (0.73 > 0.13).', 'R4'),
 ('addendum5 (+fix1 08:59 UTC)', 'r44_p1comp', 'PC1: adjusted U2 upper bound < 0; PC2: |adjusted − primary| < 0.05.', L_DC, 'Held (−0.20, upper bound −0.11); held.', 'R7'),
 ('addendum6 (28 Oct 11:50 UTC)', 'r44_redsame', 'PR1: upper bound of the reduced computation\'s share of the compression < 0.1; PR2 (descriptive prediction): overestimate of d between ensembles at least the mean within-ensemble overestimate.', L_DC, 'PR1 not held: 4% [−6, 36]; PR2 not held: 6.7 × 10⁻⁴ between, 1.3 × 10⁻³ within.', 'R5'),
 ('addendum7 (28 Oct 15:39 UTC)', 'r45_f1', '(a) C-2 lower bound > 0 with bond-direction profiles → the association is present in bond-direction fluctuations after the network-model control; (b) C-2 lower bound ≤ 0 and paired difference (whole frame − bond direction) lower bound > 0 → much of it lies outside bond directions; (c) otherwise unresolved.', L_DC, '(a): C-2 0.168 [0.071, 0.269]; association not confined to rotations about the bond.', 'R8, Fig 7'),
 ('addendum8 (29 Oct 11:15 UTC)', 'r45_mpstart', 'Descriptive (no decision rule).', L_D, 'Six ensembles × nine starts: at most 8 × 10⁻¹¹ in d.', 'Methods'),
 ('addendum9 (29 Oct 16:34 UTC; fix1 17:19 UTC)', 'r46_align', 'For A, B0 and B2: C-1 and C-2 lower bounds > 0 → holds; C-1 only → partly.', L_DC, 'Holds for A, B0 and B2.', 'R8, S3 Table'),
 ('addendum10 (30 Oct 06:13 UTC; fixes 06:46, 06:47 UTC)', 'r47_ties', 'Descriptive: largest change with average ranks and any change relative to zero.', L_D, 'At most 5 × 10⁻⁴; no interval changed its position relative to zero.', 'Methods'),
 ('addendum11 (30 Oct 07:30 UTC)', 'r48_identity_range', 'Below 15% identity: C-1 and C-2 lower bounds > 0 → holds; C-1 only → partly; comparison with the range-restriction expectation.', L_DC, 'Partly (C-1 0.167 [0.052, 0.283]; C-2 0.115 [−0.002, 0.238]); A ok, B0 not; consistent with range restriction.', 'R8, S3 Table'),
 ('addendum12 (31 Oct 13:51 UTC; fix1 14:25 UTC)', 'r49_chainbreak', 'No domain or chain of the sets not checked on 28 October has a Cα–Cα bond longer than 4.2 Å (s > 1.1).', L_SC, 'Held: none of 2,594 domains and chains (longest bond 4.149 Å).', 'Methods'),
]
CLASS = [  # section, analysis, classification, record, result
 ('R1', 'Coordinate round trip', 'pre-specified', 'round36_roundtrip 20 Oct 17:23:15 (stamp)', '20 Oct 17:26'),
 ('R1', 'Certification of 500 pairs; relative error bound', 'pre-specified', 'round36_sec6 20 Oct 17:34:46', '20 Oct 18:56'),
 ('R1', 'Computational cost', 'pre-specified', 'round41 g5 22 Oct 07:33:57', '22 Oct 08:07–08:42'),
 ('R1, S1 Fig', 'Synthetic sequences, random generator', 'pre-specified', 'round41 g6 21 Oct 08:31:54', '21 Oct 09:13'),
 ('R1, S1 Fig', 'Synthetic sequences, adversarial generators', 'added after a result', 'round41 g6_adv 21 Oct 09:25:10', '21 Oct 09:35'),
 ('R2', '12 test pairs against TM-align; factorial comparison with its d₀ sensitivity; ranking of existing measures; reversed chains (S2 Text)', 'in no record', 'before the records began (17 Oct 12:47:55)', '16 Oct'),
 ('S2 Text', 'Linear-form counterexample', 'in no record', 'before the records began (17 Oct 12:47:55)', '17 Oct 05:31'),
 ('R2', 'Chaining test (Fig 2B)', 'pre-specified', 'round17 18 Oct 00:20:07; round29 18 Oct 08:40:18', '18 Oct 00:21–00:23; 08:42'),
 ('R3', 'Calibration E0′–E4', 'pre-specified', 'round39 20 Oct 17:34:46 + amend1 21 Oct 08:27:34', '21 Oct 12:25–12:35'),
 ('R3, Table 2', 'Replicate and frame audits', 'added after a result', 'addendum 28 Oct 04:47:45', '28 Oct 05:01'),
 ('R3, Table 2', 'Earlier audits: window, survival', 'pre-specified', 'round37 18 Oct 17:08:49', '18 Oct 17:10 (survival), 17:18 (window)'),
 ('R3', 'Window difference traced to criterion and set', 'added after a result', 'round39_posthoc_window_diag 21 Oct 12:43:00', '21 Oct'),
 ('R4', '234 pairs; continuous-axis predictions', 'pre-specified', 'round36 document §1; rules in r38_b6.py, checksum in record 19 Oct 19:42:45 (after the first-session checks, which computed none of its decision quantities)', '20 Oct 01:54:50'),
 ('R4', 'Size control; W₂ against C', 'added after a result', 'addendum3 28 Oct 06:04:10; addendum5 08:53:17', '28 Oct 06:04:37; 08:55'),
 ('R4', 'κ of each cell', 'pre-specified', 'round39 E0′ 20 Oct 17:34:46', '21 Oct 12:25'),
 ('R4, S2 Table', 'Length matching', 'added after a result', 'round43 22 Oct 08:54:20', '22 Oct 09:02'),
 ('R4', 'κ² decomposition for AlphaFlow', 'added after a result', 'addendum (r44_inv split) 28 Oct 04:47:45', '28 Oct 04:48'),
 ('R5', 'E7 on 30 domains', 'pre-specified', 'round36 document §3; full text in record 19 Oct 14:56:51', '19 Oct 20:29'),
 ('R5', 'E7 on 100 domains', 'added after a result', '7g, round36_amend_C 20 Oct 06:23', '20 Oct 15:33'),
 ('R5', 'Decomposition of the compression (380 cells)', 'added after a result', 'round41 g2 21 Oct 08:31:54', '21 Oct 09:18 (CPU) and 14:44 (GPU)'),
 ('R5', 'Reduced computation on same-sequence cells', 'added after a result', 'addendum6 28 Oct 11:50:14', '28 Oct 12:29:52'),
 ('R5, Methods', 'Overestimate between sequences (at most 3.5%)', 'pre-specified', 'round36 document §0; full text in record 19 Oct 14:56:51', '19 Oct 18:55'),
 ('R6', 'Primary zone, aligned-coordinate version, secondary 1, superfamily-level 23, four domains > 320 residues removed', 'pre-specified', 'round36 7c, full text in record 18 Oct 18:58:51 (amended 19:03:50); 7d, full text in record 19 Oct 17:07:53', '19 Oct 23:03–23:06'),
 ('R6', 'Interval of secondary 1', 'added after a result', 'round43 22 Oct 08:54:20', '22 Oct 08:58'),
 ('R6', 'Folded sample counts', 'in no record', '—', '29 Oct 11:36'),
 ('R6', 'Direction, spread, profile; AlphaFlow comparison', 'pre-specified', 'round39 E5 20 Oct 17:34:46', '20 Oct 17:44'),
 ('R7', 'Primary test P1', 'pre-specified', 'round36 document §7; rules in r38_p1.py, checksum in record 20 Oct 04:58:59', '20 Oct 12:57:36'),
 ('R7', 'Power diagnostic', 'added after a result', 'round36_p1_diag 20 Oct 13:15:38 (stamp 13:15:29)', '20 Oct 13:16'),
 ('R7', 'Replication P1′ and its descriptive items', 'pre-specified', 'round36_amend_P1prime 20 Oct 13:36:03', '20 Oct 23:53'),
 ('R7', 'Composition; absolute windows', 'added after a result', 'addendum5 28 Oct 08:53:17 (+fix1 08:59:08); addendum2 05:56:13', '28 Oct 09:02; 06:13'),
 ('R8', 'First homolog study P3', 'pre-specified', 'round36 document 1b; rules in r38_p3.py, checksum in record 19 Oct 18:37:38 (snapshot copy; retained version 18:37:46)', '20 Oct 01:29:58'),
 ('R8', 'Design effect; smallest detectable effect', 'in no record', '—', '20 Oct'),
 ('R8', 'Calibration simulation (exploratory interval)', 'in no record', '—', '20 Oct 04:48–05:03'),
 ('R8', 'Confirmatory C-1, C-2; auxiliary analyses', 'pre-specified', 'round36_amend_C 20 Oct 06:23:45', '20 Oct 06:41–14:43'),
 ('R8', 'Network-model parameters', 'added after a result', 'round43 22 Oct 08:54:20', '22 Oct 09:42'),
 ('R8', 'Nonlinear and shorter-domain controls; bond directions; pair weighting', 'added after a result', 'round36_c_nl 20 Oct 15:26:31; addendum7 28 Oct 15:39:50; addendum2 05:56:13', '20 Oct 15:46; 28 Oct 17:15; 05:56'),
 ('R8', 'Other correspondences; Gly/Pro; containment; TM-align on mean points', 'added after a result', 'addendum 28 Oct 04:47:45, addendum2 05:56:13, addendum4 06:46:46', '28 Oct 06:29–07:03'),
 ('R8', 'Decomposition of the unshared part (89%)', 'in no record', '—', '21 Oct, before 09:27:41 (values written into the withdrawn round42 entry)'),
 ('R8', 'Alignment accuracy; identity range', 'added after a result', 'addendum9 29 Oct 16:34:05 (+fix1 17:19:37); addendum11 30 Oct 07:30:01', '29 Oct 17:23; 30 Oct 07:31'),
 ('Methods', 'Mean points from nine starting frames', 'added after a result', 'addendum8 29 Oct 11:15:02', '29 Oct 11:21'),
 ('Methods', 'Tie handling', 'added after a result', 'addendum10 30 Oct 06:13:48 (+fixes 06:46, 06:47)', '30 Oct 06:47'),
 ('Methods', 'Chain-break check of the calibration and homolog-study domains; package reproduction (r44_repro)', 'in no record', '—', '28 Oct (package reproduction 05:24)'),
 ('Methods', 'Chain-break check of the other sets (r49_chainbreak)', 'added after a result', 'addendum12 31 Oct 13:51:43 (+fix1 14:25:18)', '31 Oct 14:27–16:42'),
 ('Methods, S2 Text', 'Branch-and-bound check; incomplete U-statistic; upward bias of the reduced computation; alignment uncertainty; truncated chains', 'pre-specified', 'round36 document §0; full text in record 19 Oct 14:56:51', '19 Oct 16:10–19:05'),
 ('Methods, S2 Text', 'Settings and choice rule for re-measuring the bin above the bias threshold', 'added after a result', 'round36_trackB_code 19 Oct 19:42:45 (as round36_trackB_code_2026-10-20; re-stamped 19:43:30)', '19 Oct 20:39'),
 ('S2 Text', 'Diagnostics of the failed truncated-chain check', 'in no record', '—', '19 Oct 16:40–16:49'),
]
def fmt(path):
    """Column widths, wrapped text, bold frozen header."""
    from openpyxl import load_workbook
    from openpyxl.styles import Alignment, Font
    wb = load_workbook(path)
    for ws in wb.worksheets:
        for col in ws.columns:
            n = max(len(str(c.value)) if c.value is not None else 0 for c in col)
            ws.column_dimensions[col[0].column_letter].width = min(max(10, n * 0.9 + 2), 70 if ws.title not in ('Notes', 'Note') else 140)
            for c in col: c.alignment = Alignment(wrap_text=True, vertical='top')
        for c in ws[1]: c.font = Font(bold=True)
        ws.freeze_panes = 'A2'
    wb.save(path)
NAMES = {'r44_inv': 'Replicate and frame audits (r44_inv)', 'r44_corr': 'Alternative correspondences (r44_corr)', 'r44_gp': 'Glycine and proline (r44_gp)',
         'r44_pairw': 'Pair-weighted estimates (r44_pairw)', 'r44_p1win': 'Absolute trajectory windows, thermophile test (r44_p1win)',
         'r44_b6size': 'Size control across sequences (r44_b6size)', 'r44_agree': 'Agreement with TM-align on medoids (r44_agree)', 'r44_tmmean': 'TM-align on mean points (r44_tmmean)',
         'r44_w2gap': 'Gap between W₂ and Δ²_D (r44_w2gap)', 'r44_p1comp': 'Composition of thermophile and control pairs (r44_p1comp)',
         'r44_redsame': 'Reduced computation on same-sequence cells (r44_redsame)', 'r45_f1': 'Profiles restricted to bond directions (r45_f1)',
         'r45_mpstart': 'Mean points and the starting frame (r45_mpstart)', 'r46_align': 'Alignment accuracy (r46_align)', 'r47_ties': 'Tie handling (r47_ties)',
         'r48_identity_range': 'Identity below 15% (r48_identity_range)', 'r49_chainbreak': 'Chain breaks in the other sets (r49_chainbreak)', 'First session 0.0': 'Computation check 0.0: branch-and-bound against the heuristic',
         'First session 0.1 (truncated chains)': 'Computation check 0.1: truncated chains', 'First session 0.2 (upper bias of the reduced computation)': 'Computation check 0.2: upward bias of the reduced computation',
         'First session 0.3 and 0.4': 'Computation checks 0.3 and 0.4: incomplete U-statistic and flags'}
def build():
    e = entries(); a = added()
    p = pd.DataFrame(PRED, columns=['Record (UTC)', 'Analysis', 'Prediction or rule', 'Label', 'Outcome', 'Where reported'])
    p['Analysis'] = [NAMES.get(x, x) for x in p['Analysis']]
    c = pd.DataFrame(CLASS, columns=['Section', 'Analysis', 'Classification', 'Record (UTC)', 'Result (UTC)'])
    notes = pd.DataFrame({'Note': [
        'Sheet A lists every entry of the first time-stamped record (prereg_2026-10-17.json: 64 entries after two header fields, the time the record was begun, 17 October 2026 12:47:55 UTC, and a note; 101 retained Drive versions and 46 snapshot copies) with the time written inside it, the time it first appeared on Drive, later changes, in-place edits, removals, result files that already existed when it appeared, and the two flagged cases. A result file listed as existing earlier belongs to a different analysis unless the Note says otherwise. Entries of an earlier exploratory stage (17–18 October) that this paper does not use are marked; their predictions and outcomes are in the deposited record.',
        'Sheet B lists the twelve records added from 28 October 2026 onward and their five corrections; file names carry the date in Korea Standard Time (UTC+9), times inside the files are UTC.',
        'Sheet C lists the pre-specified predictions and decision rules of the analyses in this paper and of analyses specified in the records but not discussed in the text, with their outcomes. Labels: decision criterion (decides a stated reading), sanity check (expected to hold if the computation is correct), descriptive (reported without a decision), by definition (holds by construction).',
        'Sheet D classifies every analysis in the paper as pre-specified (in a record written before its result), added after a result (specified in a record after the result it follows and before it was computed) or in no record. Record times are the earliest Drive time of the entry or sub-entry (a retained version or a snapshot copy), unless marked (stamp), the time written inside it; for the records added from 28 October, the time written inside the file.',
        'Most entries hold their predictions and decision rules as text. The round36 entry (18 October 2026, 18:58:51 UTC) names a pre-specification document as the source of its sections 0–8 and holds the text of sections 7b and 7c; the record holds the full text of sections 0 and 3 from 19 October, 14:56:51. Sections 7d–7h of the document are amendments made from 19 October onward; each was entered in the record when it was made (sheet A). For sections 1, 1b and 7, the earliest evidence on Drive is the checksum of the analysis script that implements their decision rules (r38_b6.py, r38_p3.py, r38_p1.py), entered before the decision statistics were computed (for r38_p1.py, its pair-counting and mean-point steps had started up to 2.5 minutes earlier; see the Note of round36_p1_code); the scripts contain the rules and are deposited with matching checksums, and their Drive copies were last modified before those entries. For section 2 (natural insertion–deletion cell), the entry round36_r10a (20 October, 16:55:13 on Drive; stamp inside 16:43:20) holds the decision criteria of predictions 2 and 3 as text and the script checksum, before the run (19:16). The document is deposited too, but its earlier versions are dated only by local file times.',
        "In sheets C and D, R1–R8 are the eight Results sections in order (R1, Framed chains and an edit-distance metric across lengths; R2, TM-align's objective as an edit-distance metric; R3, Energy statistics for ensembles and their calibration; R4, Ensembles of different sequences on the same statistic; R5, Bond-level comparison once correspondence and rotation are fixed; R6, Errors of a temperature-conditioned generative model; R7, Thermophilic and mesophilic homologs at equal temperature; R8, Sequence identity and shared-bond fluctuation profiles). Labels such as R3′, R10a′, E1–E7 or P1 in the analysis names are the names used in the records, not section numbers.",
        'Chain-break check (Methods). A chain break would appear as a Cα–Cα bond longer than 4.2 Å (s > 1.1). On 28 October 2026 (in no record): the 336 calibration domains, every bond of the first 320 K frame of replicate 0, and the 1,496 and 949 domains of the two homolog studies, the frame-averaged bond lengths stored with each domain. On 31 October (addendum12, r49_chainbreak): the 67 model-error domains, the 494 domains of the 281 pairs of the comparison across sequences, the 34 domains of the natural insertion–deletion pairs, the 570 domains of the certified pairs and the 659 and 531 domains of the two thermophile studies, every bond of the first frame of each 320 K replicate in the mdCATH files, with the atoms each analysis used (for 311 domains of the replication, the first stored frame of the second half); the 37 force-field chains, over the residues shared with mdCATH, in the ATLAS starting structure and in mdCATH; and the 666 static chains of the first two Results sections. None had a longer bond (2,594 domains and chains on 31 October, some also checked on 28 October; longest bond 4.149 Å). Generated samples (aSAMt, AlphaFlow) were not checked: the models build them for the full domain sequence, so no residue is missing.',
        'Times of entries are the earliest Drive times at which they appear, from a retained file version or a snapshot copy (server times); times written inside entries are given where they differ.',
        'The records are included unchanged in the data deposit.']})
    with pd.ExcelWriter(f'{OUT}/S1_Table.xlsx') as xw:
        notes.to_excel(xw, index=False, sheet_name='Notes'); e.to_excel(xw, index=False, sheet_name='A Record entries'); a.to_excel(xw, index=False, sheet_name='B Added records')
        p.to_excel(xw, index=False, sheet_name='C Predictions'); c.to_excel(xw, index=False, sheet_name='D Analyses')
    fmt(f'{OUT}/S1_Table.xlsx')
    return e, a, p, c
