# Provenance

md5 of every run-time file (`code/`, `framespace_wip/experiments/`, the setup scripts and freeze records in `environment/`, `registration/`, `run_records/`),
and whether that exact checksum was recorded before or at the time of the run.

- **R** — the full md5 appears in the first time-stamped record `results/prereg_2026-10-17.json` (77 files).
  For some scripts the matching value is a later entry recorded after a correction
  (all on 2026-10-20, times written in the entries, UTC: `r38_c.py` and `r38_verify_c.py`, `round36_amend_C` `scripts_md5_v2` 06:39:35;
  `r38_r10a.py`, `round36_r10a` `amend1` `script_md5_1b` 18:08:00; `r39_r6.py`, `round36_r6` `amend1` 18:06:51; `r39_sec6.py`,
  `round36_sec6` `fix1` 19:06:04; `r38_c_nl.py`, `round36_c_nl` `amend1` 15:31:21; `r39_core.py`, `round39` `fix1` 18:25:08;
  `r38_p1n.py`, `round36_amend_P1prime` `align_restart` 16:26:08; `r38_verify_p1n.py`, the same entry's `verify_v2` 13:44:45).
- **A** — the full md5 appears in one of the round44–49 records `registration/prereg_2026-10-*_addendum*.json` (including their corrections), written before the run (17 files).
- **P** — the checksum was recorded elsewhere in the project record (checksum manifests sent with the code, a verdict output, code notes, a run log) (12 files).
- **U** — no checksum was recorded before this repository was assembled (58 files up to round43, among them 13 run queues copied from Google Drive on 2026-11-02; 17 more from round44: the run queues, the two notebook cells kept as scripts, `r44_repro.py`, `r44_p1comp2.py` and the copies of the round44 records, which were each written once on Google Drive and never modified; 3 more from round45: `r45_folded_frac.py` and the copies of the two round45 records; 10 more from round46–48 and the manuscript: the copies of the six round46–48 records, the figure and table scripts `make_figs.py`, `make_tables.py`, `make_s1.py` and the check `negtype_eta_len2.js`; 2 more from round49: the copies of its record and of its correction; 2 more on 2026-11-02: `sippl_identity_check.py`, a check of S2 Text 8.3 added after the analyses, and `e0_ht_extra.py`, the round29 E0 script, restored byte for byte from Google Drive). These are mostly
  the early CPU-phase scripts (2026-10-16 to 10-18: comparisons of TM-score forms and existing measures, ladder, C4, H-T), data-preparation
  scripts, setup scripts and freeze files. Those that were not in the local working copy were restored from
  the project storage (Google Drive) while this repository was assembled, and their sizes were checked against
  the stored files; the first record describes what these scripts computed but does not pin their bytes.

Bundles: `framespace_pkg/framespace/*.py` concatenated in sorted order has md5 `4edba596…` and
`framespace_wip/stage0/*.py` has `6703fa26…`; both prefixes are in the first record
(`tools/record_env.py` recomputes them). framespace 0.2.0 is tag `v0.2.0` of the framedchain repository (named framespace until 2026-10-26; commit `5e4fdf9`); the first record pins it
by the bundle md5, not by a commit identifier.

**Commit identifiers.** The histories of both repositories were rewritten twice, each time without changing any
file, author or date: every rewritten commit has the same file tree hash as the commit it replaces, and the tags keep
their messages and dates. On 2026-10-26 the author name was unified and trailer lines were removed from the commit
messages; on 2026-10-31 the commit messages, written in Korean until 2026-10-29, were translated into English. Each
rewrite changed every commit identifier. Identifiers in project notes map as follows (the repositories were named
framespace and framespace-paper until 2026-10-26).

| Repository | Until 2026-10-26 | 2026-10-26 to 10-31 | From 2026-10-31 | Tag | Commit |
|---|---|---|---|---|---|
| framedchain | `540c515` | `25accdb` | `81715a4` |  | framespace 0.1.0: definitions and properties of the space |
| framedchain | `ad82be1` | `137f340` | `a6a1ad1` |  | Match the confirmed definitions of the space |
| framedchain | `16be758` | `11f44e5` | `f89546f` |  | Add a PDB reader; make constants and computation parameters adjustable |
| framedchain | `31b49bf` | `62d9646` | `96c08b3` |  | Iterate the mean point to convergence |
| framedchain | `efcdd42` | `866995c` | `5e4fdf9` | `v0.2.0` | framespace 0.2.0: transport module computing the ensemble W₂ by definition |
| framedchain | `bafd810` | `676d5af` | `1162a1e` | `v0.3.0` | framespace 0.3.0: computations and estimators fixed in the paper |
| framedchain | `34b9bf4` | `694c9c2` | `2eaa00c` | `v0.3.1` | framespace 0.3.1: release preparation (MIT licence, citation metadata, English README) and remaining fixed computations |
| framedchain | `40e6b7b` | `23c7017` | `1376b7b` | `v0.3.2` | framespace 0.3.2: docstrings, comments and messages in English (no code change); SPDX licence metadata |
| framedchain | `5d5ff73` | `746631d` | `908a209` |  | framespace 0.3.3: fixes from user-side checks |
| framedchain |  | `2a943fd` | `4f5442a` | `v0.3.3` | Update commit identifiers: the v0.2.0 commit (866995c) changed with the history cleanup |
| framedchain |  | `8de7467` | `693cd79` | `v0.4.0` | framedchain 0.4.0: rename the package from framespace to framedchain |
| framedchain |  | `80c559f` | `8b23a96` | `v0.4.1` | framedchain 0.4.1: final paper title (documentation only) |
| framedchain |  | `2e219dc` | `ae5d077` |  | README: 0.4.1 gives the same values as the as-run code (framedchain-paper code/r44_repro.py, 2026-10-28; documentation only) |
| framedchain |  | `971c67b` | `3df2da1` |  | mean_point documentation: independence of the starting frame now rests on a measurement on real ensembles (6 ensembles of 153–391 bonds and 500 structures, 9 starts each, within 8e-11 after one rotation; measured with the as-run 0.2.0, 2026-10-29; documentation only) |
| framedchain-paper | `d2baa9b` | `b426e48` | `f050e4e` |  | framespace-paper 1.0.0: analysis code of the paper (as run) and reproduction documents |
| framedchain-paper | `10b6559` | `dc70dad` | `a0d0bb2` |  | STATEMENTS: state of the repositories as it is (framespace-paper still local only) |
| framedchain-paper | `60e911d` | `8eceacf` | `df9c82c` | `v1.0.0` | STATEMENTS: both repositories on GitHub (private) |
| framedchain-paper | `1072f51` | `554892c` | `cf09732` |  | Review fixes: refer to package 0.3.2, precise wording |
| framedchain-paper | `b878195` | `f177a32` | `9f35771` |  | Refer to package 0.3.3; clearer instructions for the check environment |
| framedchain-paper |  | `b566ba3` | `5d07542` |  | Update commit identifiers: map of the identifiers changed by the history cleanup |
| framedchain-paper |  | `fa8d24a` | `deaf013` |  | Paper title and package name: framedchain-paper |
| framedchain-paper |  | `6a4362e` | `c5fb5a3` |  | Final paper title and abstract |
| framedchain-paper |  | `8f2eb3e` | `5eb1571` |  | Abstract: "Likewise" before the TM-align sentence |
| framedchain-paper |  | `ad6ea70` | `3dbe570` | `v1.1.0` | round44 (2026-10-28): check scripts registered after external review, queues, record copies, documents, 1.1.0 |
| framedchain-paper |  | `600b69f` | `a3d85e8` |  | STATEMENTS: version 1.1.0 (with round44) and tags; MANIFEST.md5 |
| framedchain-paper |  | `d149c8a` | `ca361c2` |  | round44 fifth record (2026-10-28 08:53 UTC) and its correction (08:59): r44_w2p1.py (gap between W₂ and Δ²_D), r44_p1comp2.py (composition of the P1′ pairs), two record copies, documents, MANIFEST |
| framedchain-paper |  | `fa9899f` | `bf4d561` |  | round44 sixth record (2026-10-28 11:50 UTC): r44_redsame.py (reduced computation on same-sequence cells), q_r44d.sh, record copy, documents, MANIFEST |
| framedchain-paper |  | `f91c1a1` | `49c624f` | `v1.2.0` | round45 (2026-10-28 and 10-29): bond-direction profiles, starting frame of the mean point, count of folded samples; queue, two record copies, documents, 1.2.0 |
| framedchain-paper |  | `3f8f39f` | `6bb96e8` |  | docs: interactive checks copied from the session log; hardware statement corrected |
| framedchain-paper |  | `0daa4af` | `29e186b` |  | docs/INTERACTIVE_CHECKS.md: describe the source as the client-side session log (local clock) |
| framedchain-paper |  | `612de59` | `074fcda` |  | docs/INTERACTIVE_CHECKS.md: four more checks without an output file (sections 4-7) |
| framedchain-paper |  | `467ce1e` | `d70858e` |  | README: Korea Standard Time (official name) |
| framedchain-paper |  | `7ec5221` | `f47b8d0` |  | framedchain-paper 1.3.0: round46-48 checks, figure and table scripts, manuscript map |
| framedchain-paper |  | `47cd8cb` | `4dcfe47` | `v1.3.0` | docs: run dates of the figure (2026-10-30) and table (2026-10-31) scripts |
| framedchain-paper |  | `bb700a7` | `35b54b5` |  | README rewritten as a paper companion: abstract first, figures and tables to scripts, records, layout, data, environment, notes |
| framedchain-paper |  | `3485799` | `403f63f` |  | README: script naming, hardware wording |

The md5 of every file in the repository is in `MANIFEST.md5` (`md5sum -c MANIFEST.md5`).

| File | md5 | Status | Note |
|---|---|---|---|
| `code/analysis_fold208.py` | `74c43c4b10c7fea512b92cd847d7d014` | U |  |
| `code/analysis_fold208_ablation.py` | `a214900c9efa2875e9b1e1786caa645f` | U |  |
| `code/c1_core.py` | `5f00f0a17892df644e810b39109380e9` | R | full md5 in the first record |
| `code/c4_chain_extra.py` | `547b6f5e91b6b0f9d439c6929441d844` | U |  |
| `code/c4_substitution.py` | `bd280d4beba3f2ea3f9bc896b531eb43` | U |  |
| `code/chaining.py` | `41234cabeceaba3f77965574c351fb55` | U |  |
| `code/dfrechet_1674.py` | `5b6cf8af7811f887a7158eaba1f0c81c` | U |  |
| `code/e0_ht_extra.py` | `731444b3f41a1271972a8bf1d3c10305` | U | round29 E0 (spread criterion, energy distance) and H-T check; written once to Google Drive by a notebook cell on 2026-10-18 at 08:44:12 UTC, after the round29 entry (08:40:18), which names the script but not its checksum; never modified; restored byte for byte on 2026-11-02 (3,994 bytes; the code of the cell in the session log is identical) |
| `code/e0_resolution.py` | `ab9ba7700f5188a24fb72796c2794473` | U |  |
| `code/folded_analyses.py` | `d8a4a5f4643d091b70c1cbf2ca47a958` | U |  |
| `code/ht_own_alignment.py` | `140372970478f1d336a75d0924ba1dac` | U |  |
| `code/ht_theta_meso.py` | `daaa1efb5881cabc6314bc5797c067fa` | P | checksum manifest sent with the code, 2026-10-19; project code notes (atlas/code/README.md; the date of the note is not established) |
| `code/hu1_kkl_measures.py` | `cfd265881f826c627be3f3a320967944` | U |  |
| `code/hu1_rulers.py` | `ce96a02255111d3a5629fbec4243ce88` | U |  |
| `code/ladder_components.py` | `00810cb810764da3fc57302c3ffe83a6` | U |  |
| `code/ladder_ff.py` | `2a2c6a71bda54c0112758bc79420ddc0` | U |  |
| `code/ladder_ff_folded.py` | `dafaa744d5dcc6ce6b17e40ac55e774d` | U |  |
| `code/ladder_ff_single.py` | `836c8a779cd009edac6de508cca62c8c` | U |  |
| `code/ladder_half.py` | `ea16fe51c8585020b2bc94bc3dea45de` | P | checksum manifest sent with the code, 2026-10-19 |
| `code/ladder_model_error.py` | `48809b55e7b83ba446940a48fab8446d` | U |  |
| `code/ladder_model_error_folded.py` | `0c2d4e91bfcd06eca6c07152d0e9debb` | U |  |
| `code/ladder_temp_single.py` | `961d2f27c0406c9667fabb03b2d5b142` | U |  |
| `code/ladder_temp_single_v2.py` | `6420c51e5aeb3a69c51aafa1711a261d` | U |  |
| `code/make_figs.py` | `d7783270d20e7bd260710ff8b9177ab8` | U | figure script, run on Colab 2026-11-01 (Arial, text 8 pt or larger; Fig 4 redrawn 16:57 UTC with more margin; Figs 3, 4, 5 and 7 redrawn 17:45 UTC: markers per heating step in Fig 3, a key for the κ² bar in Fig 4B, negative ticks in Fig 5A, a separator in Fig 7; Figs 1, 3 and 4 redrawn 18:29 UTC and Fig 4 again 18:31 UTC: glyphs and open circles of Fig 1B and 1D drawn true to shape on axes whose x and y units differ, markers at probability 0 unclipped in Fig 3C, 'crossing' centered in its band in Fig 4B; Figs 1, 2 and 3 redrawn 22:39 UTC: arrows of the Fig 1 glyphs without end shrink, so that they scale with the term, and mathematical type for π and U in Fig 1B, D_TM in Fig 2B and n in Fig 3B; values unchanged); all figures redrawn together 2026-11-02 02:05 UTC (in the same notebook run, the values each figure function returns were saved with the manuscript files as figure_values.json): t in Fig 1D, Q in Fig 6 and ‖a‖ in S1 Fig in mathematical type, d₀(L̄) on the linear bar of Fig 2B, open diamonds in Fig 5A; Fig 4 redrawn 03:45 UTC: the κ = 1 line under the markers; all figures redrawn together 04:52 UTC with this version, the values saved again in figure_values.json and equal to those of 02:05 UTC: ticks of Fig 3C at the simulated numbers of domains (only Fig 3 changed; the other figure files are byte for byte the same); all figures redrawn together 06:07 UTC with this version in a new Colab runtime (values equal to those of 02:05 UTC): the Fig 2A row "TM-align, mean of both", the mean of TM-align's two scores (only Fig 2 changed); all figures redrawn together 17:07 UTC with this version (values unchanged): smaller markers in Fig 3A, so that the lower bar end at +28 K shows, and subscripts in ℓ₂ and ℓ₁ of S1 Fig (only Fig 3 and S1 Fig changed); reads stored results only; not specified in advance |
| `code/make_s1.py` | `f11020ff6cb866d3404df8c6a545335d` | U | S1 Table script, run on Colab 2026-10-31 and again on 2026-11-01 and 2026-11-02 (2026-11-01 16:57 UTC: percentage intervals with % on both bounds; 17:45 UTC: one timestamp style in sheet A and the label "withdrawn"; 18:29 UTC: "2,594 distinct domains and chains" in the chain-break check, and the date of the ranking in sheet D as 16–17 October; 22:45 UTC: the intersection condition and the gate in the rule of the confirmatory C-1 and C-2, with the intersection values; 2026-11-02 02:05 and 02:08 UTC: rule of the confirmatory row (equivalence and smaller-than-minimum categories; intersection values in this table only), package reproduction at 05:14 and 05:15, κ of the cells of different sequences as added after a result (round43); 2026-11-02 03:10 UTC: in sheet C, a row for C-2 in the shorter domains of the first study, ‖a‖ in the g6 controls and the full threshold rule of first-session check 0.2; in sheet D, the shorter-domain controls in the confirmatory rows, "section 0 computation checks" in the continuous-axis row and three rows classifying the coverage simulations of θ*; 03:19 UTC: records found for rows first classified as specified in no record (now "not specified in advance") — the θ* coverage simulations and the design effect reclassified as pre-specified, the other rows give where each analysis was recorded; 2026-11-02 04:07 UTC: every analysis checked against the records — the one-sided coverage simulation of θ* and the design effect pre-specified, the condition-2 simulation of the θ* test not specified in advance; every sheet D row (60 rows, new ones for analyses not covered before, among them the 499/500 description, the 450 K check of E4, where correspondence shifts occur, the single-domain resolution E0 and the negative-type computations of S3 Text 2) gives where its analysis is recorded or documented (record entry, notebook cells in docs/INTERACTIVE_CHECKS.md, or script and output file, with times); in sheet C, the two E0 predictions of round29 and round32 and the record time of the primary rule of section 7c (full text of 19 October); in sheet A, round32 marked as bearing on this paper; 2026-11-02 04:52 UTC: in sheet C, both conditions of the 450 K version of E4 with their values and the permutation reference 0.64 of round32; in sheet C, the record times of round36_roundtrip and round36_r6; in sheet D, the record time of round36_p1_code and the result times 01:29:59 of the two P3 rows; the labels "not specified in advance" and "added after a result" throughout; 2026-11-02 06:07 UTC: in the Notes, the sets the chain-break checks did not cover; in sheet C, two intervals with ties in sort order marked as not recomputed with average ranks, the remaining labels "added after a result", round36_sec6 as the record of certification prediction 2; in sheet A, the Drive time 06:51:03 of the checksums of round36_amend_C; 2026-11-02 09:40 UTC: in the Notes and in sheet C, the 2,594 chain-break checks as 1,693 distinct mdCATH domains, 37 force-field chains and 666 static chains, and no mention of the 208-structure block, which the paper does not report; 2026-11-02 17:07 UTC: in sheet C, computation checks 0.3 (a sanity check, with its median z) and 0.4 (descriptive, as in the pre-specification document) in separate rows, and the chain-break threshold as s > 1.097 (4.2 Å / 3.83 Å)); not specified in advance |
| `code/make_tables.py` | `9e13d1363bb3a942e443a2d42795952a` | U | S2 and S3 Table script, run on Colab 2026-10-31 (S3 Table again 2026-11-01 08:16 UTC with a revised note) and both again 2026-11-01 17:45 UTC (column widths and frozen headers, numeric median columns in S2, a lower bound of −0.000 in S3 now keeps its sign) and 18:29 UTC (Pairs and Superfamilies stored as numbers in the S3 xlsx; S2 note: force-field cells, and which size-matched columns are blank; csv files unchanged) and 2026-11-02 06:07 UTC (S3 Table: the four partial correlations of precision, computed under both controls, moved to the C-2 column and labeled accordingly; the S2 Table was written again with the same content); not specified in advance |
| `code/mdcath_5T_extract.py` | `a77709c6edfef1e21207070b2c959b69` | U |  |
| `code/mdcath_5T_folded.py` | `6c65c299322b4cdbf6ae447c51481fa4` | U |  |
| `code/mdcath_spread.py` | `46c7a1ea713c6ef267409b1c38067505` | U |  |
| `code/mdcath_spread_hot.py` | `2fcea3f3cde951b279fbebb97532380a` | U |  |
| `code/negtype/len2check.js` | `76e6c6c83d94ff455091a4427a5919f5` | P | checksum manifest sent with the code, 2026-10-18 |
| `code/negtype/negtype.js` | `735bb53078b445c57f9943f1c4a8ee82` | P | checksum manifest sent with the code, 2026-10-18 |
| `code/negtype/negtype_eta_len2.js` | `8fd84e58b5910171bc14e9e72161dc33` | U | η₊ sensitivity check of S3 Text 2 (2026-10-31); not specified in advance (computation) |
| `code/negtype/so3fourier.js` | `acd2130339fb0a5984ffcd5e893d77a9` | P | checksum manifest sent with the code, 2026-10-18 |
| `code/p2_gap_consensus.py` | `3984e064d12827d299af4c819381210f` | U |  |
| `code/pu_fold208.py` | `c86ac1ec102ddf89822ea71cfd8cf425` | U |  |
| `code/pu_tm_d0_sensitivity.py` | `174d65e03dfb373aa8fd6287fc2321c7` | U |  |
| `code/pu_tm_factorial.py` | `b7ccb6972ab6133561f33f3ab9030e8b` | U |  |
| `code/pu_tm_factorial_v2.py` | `d106d6368cd6e548c007da3ceb75c287` | U |  |
| `code/r37_ff.py` | `18102464b8480103e2bd83eada4055ca` | R | full md5 in the first record |
| `code/r37_protocol.py` | `9096481eaf520f9e8e89e3a395772c6f` | U |  |
| `code/r37_subst.py` | `f46284c159cfaa2fc92c940664d60823` | U |  |
| `code/r37_trunk.py` | `61148d6a8c116fe06bf4f3c3ca6f9099` | U |  |
| `code/r37_unfold.py` | `fb0b44a9b21317eb05a053cabc624bf2` | P | project code notes (atlas/code/README.md; the date of the note is not established) |
| `code/r38_7b_judge.py` | `b503a76586f04fccd111d30c6220cca7` | P | written into the 7b verdict output before reading (15:09 UTC 2026-10-19) |
| `code/r38_7b_rank.py` | `dfc78c7fbfc4c9bab672a40f2197934a` | U |  |
| `code/r38_7b_split.py` | `a9c03fcbd980031cbdc99129ca247044` | R | full md5 in the first record |
| `code/r38_af_bar.py` | `1b2e9c498836172a2aa429038a3cda2d` | U |  |
| `code/r38_alphaflow_gen.py` | `4a40fd56852ff85a5eb3d26f0cab6d0b` | R | full md5 in the first record |
| `code/r38_asamt_gen.py` | `fd3835743a2bb724869261d4e538a2b5` | R | full md5 in the first record |
| `code/r38_b6.py` | `93de102f883b1a672a2670c76c545c72` | R | full md5 in the first record |
| `code/r38_c.py` | `62de87cfd461bfbb67b5dc8eff60b703` | R | full md5 in the first record |
| `code/r38_c_fetch.py` | `dfed48968b199025ebb64810dabbc893` | R | full md5 in the first record |
| `code/r38_c_nl.py` | `1c3466f481160e124d3fb65fcea5d694` | R | full md5 in the first record |
| `code/r38_c_probe.py` | `eacce0e9f903f8efa0c562ae32f6da61` | R | full md5 in the first record |
| `code/r38_compare.py` | `f92de61603dbb677604d5907396b2f98` | R | full md5 in the first record |
| `code/r38_compare_450.py` | `5a14815eeeeca5dac41096d8d50593c1` | U |  |
| `code/r38_dff.py` | `ebd161ad81824475a7778699213f26b1` | R | full md5 in the first record |
| `code/r38_dff_D.py` | `33383f6e35b6428fa263221f5eeaca09` | R | full md5 in the first record |
| `code/r38_e7_100.py` | `458f5d96617e75baf73bcc62f06a2d9c` | R | full md5 in the first record |
| `code/r38_e7_30.py` | `8f01d240db30cb57c2d4cd0c3c1281e5` | R | full md5 in the first record |
| `code/r38_md_ref.py` | `f8b71b922256d74ba803a708a67ef756` | R | full md5 in the first record |
| `code/r38_md_ref_p90.py` | `61256ce530b3fb0ac9fb5a09d83a1fbf` | U |  |
| `code/r38_p1.py` | `58c1062c13ab8afc72c35d809eef9500` | R | full md5 in the first record |
| `code/r38_p1_diag.py` | `99e925e656799fa77d0bcfa13f134708` | R | full md5 in the first record |
| `code/r38_p1n.py` | `22511034aceb1bfc75890f99c3803308` | R | full md5 in the first record |
| `code/r38_p1n_fetch.py` | `47d9d6bb2a5d92d8fb489ce75005ce52` | R | full md5 in the first record |
| `code/r38_p2_D.py` | `7bf20103014e50185aad942ce07bcef4` | R | full md5 in the first record |
| `code/r38_p2_D_af.py` | `b1c0f68e45eb20d01285e7ab5098038b` | R | full md5 in the first record |
| `code/r38_p2_aggregate.py` | `bd8657fdba6eabe318ac2f428c5a9222` | R | full md5 in the first record |
| `code/r38_p2_aggregate_D.py` | `de6a60a2f48cb24f29cac5d41b938ab5` | R | full md5 in the first record |
| `code/r38_p2_aggregate_D2.py` | `2c20fb3c21d23941964fa7f105f839aa` | R | full md5 in the first record |
| `code/r38_p3.py` | `a37a96530e3d5a1bd42b2691d137278a` | R | full md5 in the first record |
| `code/r38_p3_explore.py` | `4f4f0b60522677212ef274015dd276b5` | U |  |
| `code/r38_p3_explore2.py` | `ff85983485002a5eb2986e0c61d62e94` | U |  |
| `code/r38_p3_sens.py` | `1c7697146435fcb4a4dbd5bc394c5e7a` | R | full md5 in the first record |
| `code/r38_prep.py` | `380a80cd539a6a0d0dc395e9768d92d2` | U |  |
| `code/r38_r10a.py` | `76294c0b61fea343ef0c4ea31287264d` | R | full md5 in the first record |
| `code/r38_s0.py` | `504a48c0eeb52f6bbb9ca43be55938ee` | R | full md5 in the first record |
| `code/r38_s0_bias2.py` | `61229439eea8ba434e026325775aa8fd` | R | full md5 in the first record |
| `code/r38_static.py` | `9deb1c9f96de257aeb118d13aaf247f8` | R | full md5 in the first record |
| `code/r38_verify_b.py` | `beb68f06ad965a3c9a8abd6c4567db95` | R | full md5 in the first record |
| `code/r38_verify_c.py` | `784aef82b0637c344534f27ef2649e39` | R | full md5 in the first record |
| `code/r38_verify_p1.py` | `22f4e1d568e6455c82d72d98e4cc0539` | R | full md5 in the first record |
| `code/r38_verify_p1n.py` | `1b86e96f708079c4c2916a2d3e14667c` | R | full md5 in the first record |
| `code/r39_D_stem.py` | `2443e391c25a57af588ab6631208e784` | R | full md5 in the first record |
| `code/r39_D_unfold.py` | `26c22d770668b0a8f905fa1551f05c2c` | R | full md5 in the first record |
| `code/r39_core.py` | `04e57acc63848d37fb57cae3a6bb02d0` | R | full md5 in the first record |
| `code/r39_e5prep.py` | `aee877c8bd8a5d14717cf6262b84928b` | R | full md5 in the first record |
| `code/r39_e6.py` | `6729de931b792dd6d031ca57992e9402` | R | full md5 in the first record |
| `code/r39_r10a_posthoc.py` | `5f20e8591412607000320be0edf7b561` | R | full md5 in the first record |
| `code/r39_r6.py` | `a12e10a5095640c077badd83f93da713` | R | full md5 in the first record |
| `code/r39_roundtrip.py` | `a4984f12725e509edbd908fa7ad4fff4` | R | full md5 in the first record |
| `code/r39_sec6.py` | `1b11b675bfd25b2e34dafebc217b5618` | R | full md5 in the first record |
| `code/r41_compute.py` | `ac8fc3ce936e13216552da71df8bcb8d` | R | full md5 in the first record |
| `code/r41_compute_prep.py` | `add2ad62ab44ea5fc513a766a00bd7a2` | R | full md5 in the first record |
| `code/r41_dcomp.py` | `eff667970195cb45117f89b6fc9fbca5` | R | full md5 in the first record |
| `code/r41_dff50.py` | `88360ac232cad8a36a207a741f3bcdb9` | R | full md5 in the first record |
| `code/r41_gapinfo.py` | `b823de6670721190025b506f3a79b475` | R | full md5 in the first record |
| `code/r41_synth.py` | `982f6f07ac0f92d342e95b98f04cd70a` | R | full md5 in the first record |
| `code/r41_synth_adv.py` | `96e6bc9e16cfbcc80512baa2d020fd0c` | R | full md5 in the first record |
| `code/r43_sens.py` | `1643ce8b2bc0133a4a7cd37b4a16de60` | R | full md5 in the first record |
| `code/r44_b6size_cell26.py` | `41fb4e7319e297133e2dcc0e16ce1a33` | U | notebook cell 26 (2026-10-28 06:04:37 UTC) verbatim below a header added for this repository |
| `code/r44_corr.py` | `4c23a80cca4644bd52a03ffc4f3b7d19` | A | addendum 04:47:45 UTC |
| `code/r44_corr2.py` | `972af7945c56590d4d8ef15f964a5ce5` | P | md5 prefix written to the run log `results/logs/r44.log` at the start of the queue, 2026-10-28 06:50:45 UTC (record addendum4 names the script) |
| `code/r44_corr_v2.py` | `f364e242d17ae36075c5c761a0d76fba` | A | addendum2 05:56:13 UTC |
| `code/r44_inv.py` | `b4562167ff541206c287f9d43282c79a` | A | addendum 04:47:45 UTC |
| `code/r44_p1comp2.py` | `06b23792e1e922755be61e70678f5c25` | U | replaces part p1comp of `r44_w2p1.py` as described in `_addendum5_fix1.json`; md5 printed at upload before the run (2026-10-28 09:02 UTC) |
| `code/r44_p1win.py` | `90b241be249fe0c2f94d1d0f8cf4a8d4` | A | addendum2 05:56:13 UTC |
| `code/r44_pairw.py` | `208d60feb7ef118952d0c245f279497c` | A | addendum2 05:56:13 UTC |
| `code/r44_redsame.py` | `2a6a20cf5687a08611518b0454174a3f` | P | md5 prefix written to the run log `results/logs/r44.log` at the start of the queue, 2026-10-28 11:55:58 UTC (record addendum6 names the script); two earlier starts stopped before computing (missing package; a variable-name clash fixed in this version) |
| `code/r44_repro.py` | `a8bd5a23aa5f896d1e267dac160fedc0` | U | implementation check, run on 2026-10-28 |
| `code/r44_seqid_cell30.py` | `baba0c25848e601ce7207da39d4fa636` | U | notebook cell 30 verbatim below a header added for this repository; not specified in advance |
| `code/r44_w2p1.py` | `70d6d7f49b3847a43f947192a4506f58` | P | md5 prefix written into the correction record `registration/prereg_2026-10-28_addendum5_fix1.json` (08:59:08 UTC, after part w2 had run at 08:55 UTC); part w2 ran as recorded, part p1comp stopped before producing values |
| `code/r45_f1.py` | `55889bfa3aa916f3cabdf6d5517434b6` | A | addendum7 2026-10-28 15:39:50 UTC (`script_md5`) |
| `code/r45_folded_frac.py` | `8420f85172445aaada47480d2e7a89e8` | U | notebook run of `/content/folded_frac.py` (2026-10-29 11:31–11:37 UTC) verbatim below a header added for this repository; not specified in advance |
| `code/r45_mpstart.py` | `642e12e591c07a3f9d907a96166f05a1` | A | addendum8 2026-10-29 11:15:02 UTC (`script_md5`) |
| `code/r46_align.py` | `a2c0156b3068202f3acf1dba8ca69035` | A | addendum9 2026-10-29 16:34:05 UTC (`script_md5`) |
| `code/r46_align_fix1.py` | `6353ddb70d0187e0462d41e1e8a6a44a` | A | addendum9_fix1 2026-10-29 17:19:37 UTC (`script_md5`) |
| `code/r47_combine_fix1.py` | `b18f1d038685dc0a107824a6c7f1a58f` | A | addendum10_fix1 2026-10-30 06:46:19 UTC (`script_md5`); stopped with a syntax error, no output |
| `code/r47_combine_fix2.py` | `10ab7f1ecb2fe3ee2cec531e6794fd16` | A | addendum10_fix2 2026-10-30 06:47:24 UTC (`script_md5`) |
| `code/r47_ties.py` | `62736bf92a18c09a4b5fede7ad2ae023` | A | addendum10 2026-10-30 06:13:48 UTC (`script_md5`) |
| `code/r48_identity_range.py` | `32f24903bd78c607ba65dc41dae98d66` | A | addendum11 2026-10-30 07:30:01 UTC (`script_md5`) |
| `code/r49_chainbreak.py` | `4e5e398d0d59be9f759af402c70e0e1a` | A | addendum12 2026-10-31 13:51:43 UTC (`script_md5`); its run stopped when the Colab VM was replaced (about 14:13 UTC), with no output |
| `code/r49_chainbreak_fix1.py` | `363b9141170533127c235412d1382a1d` | A | addendum12_fix1 2026-10-31 14:25:18 UTC (`script_md5`); the same computation, resumable |
| `code/sippl_identity_check.py` | `2134b09e8f0bcfdb1cab980aeb1c07c9` | U | check of the comparison in S2 Text 8.3 (the linear-form counterexample against the identity correspondence), added on 2026-11-02 and run locally (output in its header); not specified in advance |
| `environment/af_bar_readmeenv.sh` | `b4af12134c383fbc056f2ad6ab48ffd0` | U |  |
| `environment/afenv_freeze.txt` | `0087a3eadc056c67bb37d64430650f52` | U |  |
| `environment/samenv_freeze.txt` | `046b712106ab93314518c9622d2bed40` | U |  |
| `environment/setup_afenv3.sh` | `669392c836ac336f4eafa4afe60d1c8a` | U |  |
| `environment/setup_samenv3.sh` | `d2f64ac3a09bf953ef8fd2fcf10f07fa` | U |  |
| `framespace_wip/experiments/sippl_linear_counterexample.py` | `321f399b325a2fff79d31425321621ff` | P | checksum manifest sent with the code, 2026-10-17 |
| `registration/prereg_2026-10-28_addendum.json` | `b3b5c040fd075bab96a4036bc3b8b8c9` | U | byte copy of the Drive file (created 2026-10-28 04:47:45 UTC, never modified) |
| `registration/prereg_2026-10-28_addendum2.json` | `79258a06e372c0b3fecf68047e53e8b1` | U | byte copy of the Drive file (created 05:56:13 UTC, never modified) |
| `registration/prereg_2026-10-28_addendum3.json` | `83230432bc1c7a979ce276f1ee24c90f` | U | byte copy of the Drive file (created 06:04:11 UTC, never modified) |
| `registration/prereg_2026-10-28_addendum4.json` | `2b53341c96605dc94f04c4f03abf4eab` | U | byte copy of the Drive file (created 06:46:46 UTC, never modified) |
| `registration/prereg_2026-10-28_addendum5.json` | `c0ca817cd5378cfc2934535041401677` | U | byte copy of the Drive file (created 08:53:17 UTC, never modified) |
| `registration/prereg_2026-10-28_addendum5_fix1.json` | `126fbb4d4f7b6d59b45e51cb77e79c9a` | U | byte copy of the Drive file (created 08:59:08 UTC, never modified) |
| `registration/prereg_2026-10-28_addendum6.json` | `9b4ad04cc1fd1dd74e946df7cbfbd9e8` | U | byte copy of the Drive file (created 11:50:14 UTC, never modified) |
| `registration/prereg_2026-10-29_addendum7.json` | `5c623bf34007f527c4b5c15441249908` | U | byte copy of the Drive file (written 2026-10-28 15:39:50 UTC, never modified; the file-name date is KST) |
| `registration/prereg_2026-10-29_addendum8.json` | `3d28986f75ccd406e02f4006336a673e` | U | byte copy of the Drive file (written 2026-10-29 11:15:02 UTC, never modified) |
| `registration/prereg_2026-10-30_addendum10.json` | `e2aba4016277859b350287d2696065ca` | U | byte copy of the Drive file (written 2026-10-30 06:13:48 UTC, never modified) |
| `registration/prereg_2026-10-30_addendum10_fix1.json` | `90d7f433598593c77aaa6410fa3cd4a1` | U | byte copy of the Drive file (written 2026-10-30 06:46:19 UTC, never modified) |
| `registration/prereg_2026-10-30_addendum10_fix2.json` | `393880a87a6dc028e215ecdea559c5b9` | U | byte copy of the Drive file (written 2026-10-30 06:47:24 UTC, never modified) |
| `registration/prereg_2026-10-30_addendum11.json` | `947aa2fef41764e6e6b377ee93cb4b3d` | U | byte copy of the Drive file (written 2026-10-30 07:30:01 UTC, never modified) |
| `registration/prereg_2026-10-30_addendum9.json` | `c4853aaed96544685a940df29eee4478` | U | byte copy of the Drive file (written 2026-10-29 16:34:05 UTC, never modified) |
| `registration/prereg_2026-10-30_addendum9_fix1.json` | `b61187dacccbaab2fa46ca3aff4ad940` | U | byte copy of the Drive file (written 2026-10-29 17:19:37 UTC, never modified) |
| `registration/prereg_2026-10-31_addendum12.json` | `b2c8ebcb40b986db06f640ed1098e661` | U | byte copy of the Drive file (written 2026-10-31 13:51:43 UTC, never modified) |
| `registration/prereg_2026-10-31_addendum12_fix1.json` | `f4fae98ae64fea89ab6fca63578927aa` | U | byte copy of the Drive file (written 2026-10-31 14:25:18 UTC, never modified) |
| `registration/prereg_round34.json` | `0453e32ea141d0c8d6d1258542e3ebd8` | P | checksum manifest sent with the code, 2026-10-19 |
| `registration/prereg_round37.json` | `3692c5e053c1a3fd87a22df2eefa34e7` | U |  |
| `run_records/queues/q_7b.sh` | `0f5807d1566cbe077b816533639da621` | R | full md5 in the first record (`round36_partial_freeze_7b7c`, `code_md5`, 2026-10-19 14:56:51 UTC); copied byte for byte from Google Drive on 2026-11-02 |
| `run_records/queues/q_7b_rank.sh` | `faf7253f4c6b01ae32b21a1c1188d2c4` | U | copied byte for byte from Google Drive (`code/q_7b_rank.sh`, file time 2026-10-19 17:10:38 UTC) on 2026-11-02 |
| `run_records/queues/q_af_par.sh` | `3d73dbc6f3cadf8bec167deef819b5f9` | U | copied byte for byte from Google Drive (`code/q_af_par.sh`, file time 2026-10-19 09:24:51 UTC) on 2026-11-02 |
| `run_records/queues/q_agg.restore_20261020.sh` | `374d8886a7fe5bd6c9a703995daf6191` | U | copied byte for byte from Google Drive (`code/q_agg.restore_20261020.sh`, file time 2026-10-20 06:02:24 UTC) on 2026-11-02; the queue restarted after the Colab VM was replaced on 2026-10-20 (it waits for the restored environment), with the same commands as `q_agg.sh` |
| `run_records/queues/q_agg.sh` | `152e002a5b8d2fba86cd585980bbb1d2` | R | full md5 in the first record (`round36_partial_freeze_7b7c`, `code_md5`, 2026-10-19 14:56:51 UTC); copied byte for byte from Google Drive on 2026-11-02 |
| `run_records/queues/q_b6_cpu.sh` | `ab6ed598d287f688b1205816f6ec7431` | R | full md5 in the first record |
| `run_records/queues/q_c_cpu.sh` | `a83e8fd847115e5c4338d26eed3ee592` | R | full md5 in the first record |
| `run_records/queues/q_compare_cpu.restore_20261020.sh` | `11a07afdcbf714249fd335a0d2f735a9` | U | copied byte for byte from Google Drive (`code/q_compare_cpu.restore_20261020.sh`, file time 2026-10-20 06:02:24 UTC) on 2026-11-02; the queue restarted after the Colab VM was replaced on 2026-10-20 (it waits for the restored environment), with the same commands as `q_compare_cpu.sh` |
| `run_records/queues/q_compare_cpu.sh` | `a8586d9d1aa7c42cf9c2e664ab482c4d` | R | full md5 in the first record (`round36_partial_freeze_7b7c`, `code_md5`, 2026-10-19 14:56:51 UTC); copied byte for byte from Google Drive on 2026-11-02 |
| `run_records/queues/q_ct.sh` | `2b06e77aa4daf92a60b673b2459681bf` | U | copied byte for byte from Google Drive (`code/q_ct.sh`, file time 2026-10-22 07:41:11 UTC) on 2026-11-02 |
| `run_records/queues/q_full_gpu.sh` | `cfc6fa2417882c0d844ac3192fe79099` | U | copied byte for byte from Google Drive (`code/q_full_gpu.sh`, file time 2026-10-19 09:41:49 UTC) on 2026-11-02 |
| `run_records/queues/q_full_gpu10.sh` | `5dfe38f138ea1f8f221cc3e568155af4` | R | full md5 in the first record |
| `run_records/queues/q_full_gpu11.sh` | `34167bc64c353fad89e1dcbbdac10c62` | R | full md5 in the first record (`round36_scope_2026-10-20`, `queues`, 16:39:16 UTC); copied byte for byte from Google Drive on 2026-11-02 |
| `run_records/queues/q_full_gpu12.sh` | `08d02727207e26b7ad4e89e7953aac69` | R | full md5 in the first record (`round36_r10a_2026-10-20`, queue, 16:55:02 UTC); until 2026-11-02 the copy here, retyped from the session log, had one more final newline |
| `run_records/queues/q_full_gpu13.sh` | `ccc8247223c973479a9eab26bf93fe9d` | R | full md5 in the first record |
| `run_records/queues/q_full_gpu14.sh` | `84367c203520890f48c111d4150439a7` | R | full md5 in the first record |
| `run_records/queues/q_full_gpu2.sh` | `ca26e968ad230332050fe8f583aef190` | R | full md5 in the first record (`round36_partial_freeze_7b7c`, `code_md5`, 2026-10-19 14:56:51 UTC); copied byte for byte from Google Drive on 2026-11-02 |
| `run_records/queues/q_full_gpu3.sh` | `83b70fb51cfcaeb5e4efe34b28022759` | U | copied byte for byte from Google Drive (`code/q_full_gpu3.sh`, file time 2026-10-19 16:58:29 UTC) on 2026-11-02 |
| `run_records/queues/q_full_gpu4.sh` | `369abb4d5d22b0e2fa58d544dc0966b1` | R | full md5 in the first record (`round36_amend_p2_D_2026-10-19`, `scripts_md5`, 17:11:44 UTC); copied byte for byte from Google Drive on 2026-11-02 |
| `run_records/queues/q_full_gpu5.sh` | `f6b552ff5a1ea7a91cb0040ff8f4f1aa` | R | full md5 in the first record (`round36_p3_spec_2026-10-19`, `queue_md5`, 18:37:38 UTC); copied byte for byte from Google Drive on 2026-11-02 |
| `run_records/queues/q_full_gpu6.sh` | `d1c6079e97024d7a24eaf8ee8ba2f8e4` | R | full md5 in the first record |
| `run_records/queues/q_full_gpu7.sh` | `df44ecd21aa844400c0546912199cc4a` | U | copied byte for byte from Google Drive (`code/q_full_gpu7.restore_20261020.sh`, file time 2026-10-20 06:02:24 UTC) on 2026-11-02; the queue restarted after the Colab VM was replaced on 2026-10-20 (it waits for the restored environment) |
| `run_records/queues/q_full_gpu8.sh` | `3247137d47931d6a21d3dec884bec206` | R | full md5 in the first record |
| `run_records/queues/q_full_gpu9.sh` | `f2cd9db4e000eaa36bef403178ca751f` | R | full md5 in the first record |
| `run_records/queues/q_p1_cpu.restore_20261020.sh` | `12059d0a1d3df48097ed7b4f928cf226` | U | copied byte for byte from Google Drive (`code/q_p1_cpu.restore_20261020.sh`, file time 2026-10-20 06:02:24 UTC) on 2026-11-02; the queue restarted after the Colab VM was replaced on 2026-10-20 (it waits for the restored environment), with 2 workers instead of 3 |
| `run_records/queues/q_p1_cpu.sh` | `98dafc7db4669fa0c867363d569e51a6` | R | full md5 in the first record |
| `run_records/queues/q_p1n.sh` | `324a38025b2cd24d84f30a23f1190cbc` | R | full md5 in the first record (`round36_amend_P1prime_2026-10-20`, `scripts_md5` (as `q_p1n.restore_20261020.sh`), 13:35:54 UTC); copied byte for byte from Google Drive on 2026-11-02 |
| `run_records/queues/q_p1n2.sh` | `445c532a66b376b23c662d929233a9d1` | R | full md5 in the first record (`round36_amend_P1prime_2026-10-20`, `align_restart`, 16:26:08 UTC); copied byte for byte from Google Drive on 2026-11-02 |
| `run_records/queues/q_p3_cpu.sh` | `c3eaa91d95bcac070020ab0dd4411271` | R | full md5 in the first record (`round36_p3_spec_2026-10-19`, `queue_md5`, 18:37:38 UTC); copied byte for byte from Google Drive on 2026-11-02 |
| `run_records/queues/q_p3s_cpu.sh` | `dc9952e7a634a33e7d9489cbb5644e0c` | R | full md5 in the first record |
| `run_records/queues/q_r39_cpu1.sh` | `138ad916e3faa270352f6062a6559f40` | R | full md5 in the first record |
| `run_records/queues/q_r39_cpu2.sh` | `a5993779ffc78a505723288a2770853c` | R | full md5 in the first record |
| `run_records/queues/q_r39_e6.sh` | `9839e6335dd8bc5e8f32ef4156394e23` | U | copied byte for byte from Google Drive (`code/q_r39_e6.sh`, file time 2026-10-21 09:36:31 UTC) on 2026-11-02 |
| `run_records/queues/q_r39_gpu.sh` | `9a0e52c775f86001fab404835bf97dca` | R | full md5 in the first record |
| `run_records/queues/q_r39_gpu2.sh` | `b8f29001152d1d10b9719187cfe45c9a` | R | full md5 in the first record (`round36_scope_2026-10-20`, `queues`, 16:39:16 UTC); copied byte for byte from Google Drive on 2026-11-02 |
| `run_records/queues/q_r41_cpu1.sh` | `587094db31289ed9147715632dec2321` | U | copied byte for byte from Google Drive (`code/q_r41_cpu1.sh`, file time 2026-10-21 08:36:31 UTC) on 2026-11-02 |
| `run_records/queues/q_r41_cpu2.sh` | `39691ef944c26a532a5b2844d6010b6a` | U | copied byte for byte from Google Drive (`code/q_r41_cpu2.sh`, file time 2026-10-21 08:42:46 UTC) on 2026-11-02 |
| `run_records/queues/q_r41_gpu.sh` | `1950456826e8af5e330a5a0e48345832` | U | copied byte for byte from Google Drive (`code/q_r41_gpu.sh`, file time 2026-10-21 08:46:14 UTC) on 2026-11-02 |
| `run_records/queues/q_r44_cpu.sh` | `7c43fc43d4263cdfb24a65da8b02419e` | U |  |
| `run_records/queues/q_r44_gpu.sh` | `c535049adaa9335d4e4422be543dc394` | U |  |
| `run_records/queues/q_r44b_corr.sh` | `12e058930ddb2843abdc0b856fcb9ece` | U |  |
| `run_records/queues/q_r44b_misc.sh` | `7112c17b8c2e46a5a17b29b8440f5cbe` | U |  |
| `run_records/queues/q_r44c.sh` | `0a582e7eab7640bb6b9e67fad46543df` | U |  |
| `run_records/queues/q_r44d.sh` | `c36712c2119c459dcfaa25a13eb0ced5` | U |  |
| `run_records/queues/q_r45.sh` | `3b1de9082ceff5a5b51ce9f6cd2071da` | A | addendum7 (`queue_md5`) |
| `run_records/queues/q_r47.py` | `d2b027c166d970df52a0dcdbec187a7c` | A | addendum10 (`runner_md5`) |
| `run_records/queues/r39_resume_20261021.sh` | `a32e0766f12b37572dbc215c9dbe295c` | U |  |
