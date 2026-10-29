# Provenance

md5 of every run-time file (`code/`, `framespace_wip/experiments/`, the setup scripts and freeze records in `environment/`, `registration/`, `run_records/`),
and whether that exact checksum was recorded before or at the time of the run.

- **R** — the full md5 appears in the registration record `results/prereg_2026-10-17.json` (65 files).
  For some scripts the matching value is a later entry registered after a correction
  (`r38_c.py`, `r38_verify_c.py`, `r38_r10a.py`, `r39_r6.py`, `r39_sec6.py`); the record gives its time.
- **A** — the full md5 appears in one of the round44 or round45 records `registration/prereg_2026-10-2*_addendum*.json`, written before the run (8 files).
- **P** — the checksum was recorded elsewhere in the project record (checksum manifests sent with the code, a verdict output, code notes, a run log) (11 files).
- **U** — no checksum was recorded before this repository was assembled (46 files up to round43; 18 more from round44: the run queues, the two notebook cells kept as scripts, `r44_repro.py`, `r44_p1comp2.py`, `r44_redsame.py` and the copies of the round44 records, which were each written once on Google Drive and never modified; 3 more from round45: `r45_folded_frac.py` and the copies of the two round45 records). These are mostly
  the early CPU-phase scripts (2026-10-16 to 10-18: ruler comparisons, ladder, C4, H-T), data-preparation
  scripts, setup scripts and freeze files. Those that were not in the local working copy were restored from
  the project storage (Google Drive) while this repository was assembled, and their sizes were checked against
  the stored files; the registration record describes what these scripts computed but does not pin their bytes.

Bundles: `framespace_pkg/framespace/*.py` concatenated in sorted order has md5 `4edba596…` and
`framespace_wip/stage0/*.py` has `6703fa26…`; both prefixes are in the registration record
(`tools/record_env.py` recomputes them). framespace 0.2.0 is tag `v0.2.0` of the framedchain repository (named framespace until 2026-10-26; commit `866995c`); the registration record pins it
by the bundle md5, not by a commit identifier.

**Commit identifiers.** On 2026-10-26 the author name in the histories of both repositories was unified and trailer
lines were removed from the commit messages. This changed every commit identifier but no file or date: each
rewritten commit has the same file tree hash as the commit it replaces. Identifiers in project notes written before that date map as follows (the repositories are now named framedchain
and framedchain-paper).

| Repository | Before | Now | Tag | Commit |
|---|---|---|---|---|
| framespace | `540c515` | `25accdb` | | framespace 0.1.0 |
| framespace | `ad82be1` | `137f340` | | definitions aligned with the space document |
| framespace | `16be758` | `11f44e5` | | PDB reader, adjustable constants |
| framespace | `31b49bf` | `62d9646` | | mean point iterated to convergence |
| framespace | `efcdd42` | `866995c` | `v0.2.0` | framespace 0.2.0, the version used in the runs (file tree `6db30c0`) |
| framespace | `bafd810` | `676d5af` | `v0.3.0` | framespace 0.3.0 |
| framespace | `34b9bf4` | `694c9c2` | `v0.3.1` | framespace 0.3.1 |
| framespace | `40e6b7b` | `23c7017` | `v0.3.2` | framespace 0.3.2 |
| framespace | `5d5ff73` | `746631d` | | framespace 0.3.3 (tag `v0.3.3` is on the next commit, `2a943fd`, which updates the cited identifiers) |
| framespace-paper | `d2baa9b` | `b426e48` | | framespace-paper 1.0.0 |
| framespace-paper | `10b6559` | `dc70dad` | | |
| framespace-paper | `60e911d` | `8eceacf` | `v1.0.0` | first push |
| framespace-paper | `1072f51` | `554892c` | | |
| framespace-paper | `b878195` | `f177a32` | | |

The md5 of every file in the repository is in `MANIFEST.md5` (`md5sum -c MANIFEST.md5`).

| File | md5 | Status | Note |
|---|---|---|---|
| `code/analysis_fold208.py` | `74c43c4b10c7fea512b92cd847d7d014` | U |  |
| `code/analysis_fold208_ablation.py` | `a214900c9efa2875e9b1e1786caa645f` | U |  |
| `code/c1_core.py` | `5f00f0a17892df644e810b39109380e9` | R | full md5 in the registration record |
| `code/c4_chain_extra.py` | `547b6f5e91b6b0f9d439c6929441d844` | U |  |
| `code/c4_substitution.py` | `bd280d4beba3f2ea3f9bc896b531eb43` | U |  |
| `code/chaining.py` | `41234cabeceaba3f77965574c351fb55` | U |  |
| `code/dfrechet_1674.py` | `5b6cf8af7811f887a7158eaba1f0c81c` | U |  |
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
| `code/mdcath_5T_extract.py` | `a77709c6edfef1e21207070b2c959b69` | U |  |
| `code/mdcath_5T_folded.py` | `6c65c299322b4cdbf6ae447c51481fa4` | U |  |
| `code/mdcath_spread.py` | `46c7a1ea713c6ef267409b1c38067505` | U |  |
| `code/mdcath_spread_hot.py` | `2fcea3f3cde951b279fbebb97532380a` | U |  |
| `code/negtype/len2check.js` | `76e6c6c83d94ff455091a4427a5919f5` | P | checksum manifest sent with the code, 2026-10-18 |
| `code/negtype/negtype.js` | `735bb53078b445c57f9943f1c4a8ee82` | P | checksum manifest sent with the code, 2026-10-18 |
| `code/negtype/so3fourier.js` | `acd2130339fb0a5984ffcd5e893d77a9` | P | checksum manifest sent with the code, 2026-10-18 |
| `code/p2_gap_consensus.py` | `3984e064d12827d299af4c819381210f` | U |  |
| `code/pu_fold208.py` | `c86ac1ec102ddf89822ea71cfd8cf425` | U |  |
| `code/pu_tm_d0_sensitivity.py` | `174d65e03dfb373aa8fd6287fc2321c7` | U |  |
| `code/pu_tm_factorial.py` | `b7ccb6972ab6133561f33f3ab9030e8b` | U |  |
| `code/pu_tm_factorial_v2.py` | `d106d6368cd6e548c007da3ceb75c287` | U |  |
| `code/r37_ff.py` | `18102464b8480103e2bd83eada4055ca` | R | full md5 in the registration record |
| `code/r37_protocol.py` | `9096481eaf520f9e8e89e3a395772c6f` | U |  |
| `code/r37_subst.py` | `f46284c159cfaa2fc92c940664d60823` | U |  |
| `code/r37_trunk.py` | `61148d6a8c116fe06bf4f3c3ca6f9099` | U |  |
| `code/r37_unfold.py` | `fb0b44a9b21317eb05a053cabc624bf2` | P | project code notes (atlas/code/README.md; the date of the note is not established) |
| `code/r38_7b_judge.py` | `b503a76586f04fccd111d30c6220cca7` | P | written into the 7b verdict output before reading (15:09 UTC 2026-10-19) |
| `code/r38_7b_rank.py` | `dfc78c7fbfc4c9bab672a40f2197934a` | U |  |
| `code/r38_7b_split.py` | `a9c03fcbd980031cbdc99129ca247044` | R | full md5 in the registration record |
| `code/r38_af_bar.py` | `1b2e9c498836172a2aa429038a3cda2d` | U |  |
| `code/r38_alphaflow_gen.py` | `4a40fd56852ff85a5eb3d26f0cab6d0b` | R | full md5 in the registration record |
| `code/r38_asamt_gen.py` | `fd3835743a2bb724869261d4e538a2b5` | R | full md5 in the registration record |
| `code/r38_b6.py` | `93de102f883b1a672a2670c76c545c72` | R | full md5 in the registration record |
| `code/r38_c.py` | `62de87cfd461bfbb67b5dc8eff60b703` | R | full md5 in the registration record |
| `code/r38_c_fetch.py` | `dfed48968b199025ebb64810dabbc893` | R | full md5 in the registration record |
| `code/r38_c_nl.py` | `1c3466f481160e124d3fb65fcea5d694` | R | full md5 in the registration record |
| `code/r38_c_probe.py` | `eacce0e9f903f8efa0c562ae32f6da61` | R | full md5 in the registration record |
| `code/r38_compare.py` | `f92de61603dbb677604d5907396b2f98` | R | full md5 in the registration record |
| `code/r38_compare_450.py` | `5a14815eeeeca5dac41096d8d50593c1` | U |  |
| `code/r38_dff.py` | `ebd161ad81824475a7778699213f26b1` | R | full md5 in the registration record |
| `code/r38_dff_D.py` | `33383f6e35b6428fa263221f5eeaca09` | R | full md5 in the registration record |
| `code/r38_e7_100.py` | `458f5d96617e75baf73bcc62f06a2d9c` | R | full md5 in the registration record |
| `code/r38_e7_30.py` | `8f01d240db30cb57c2d4cd0c3c1281e5` | R | full md5 in the registration record |
| `code/r38_md_ref.py` | `f8b71b922256d74ba803a708a67ef756` | R | full md5 in the registration record |
| `code/r38_md_ref_p90.py` | `61256ce530b3fb0ac9fb5a09d83a1fbf` | U |  |
| `code/r38_p1.py` | `58c1062c13ab8afc72c35d809eef9500` | R | full md5 in the registration record |
| `code/r38_p1_diag.py` | `99e925e656799fa77d0bcfa13f134708` | R | full md5 in the registration record |
| `code/r38_p1n.py` | `22511034aceb1bfc75890f99c3803308` | R | full md5 in the registration record |
| `code/r38_p1n_fetch.py` | `47d9d6bb2a5d92d8fb489ce75005ce52` | R | full md5 in the registration record |
| `code/r38_p2_D.py` | `7bf20103014e50185aad942ce07bcef4` | R | full md5 in the registration record |
| `code/r38_p2_D_af.py` | `b1c0f68e45eb20d01285e7ab5098038b` | R | full md5 in the registration record |
| `code/r38_p2_aggregate.py` | `bd8657fdba6eabe318ac2f428c5a9222` | R | full md5 in the registration record |
| `code/r38_p2_aggregate_D.py` | `de6a60a2f48cb24f29cac5d41b938ab5` | R | full md5 in the registration record |
| `code/r38_p2_aggregate_D2.py` | `2c20fb3c21d23941964fa7f105f839aa` | R | full md5 in the registration record |
| `code/r38_p3.py` | `a37a96530e3d5a1bd42b2691d137278a` | R | full md5 in the registration record |
| `code/r38_p3_explore.py` | `4f4f0b60522677212ef274015dd276b5` | U |  |
| `code/r38_p3_explore2.py` | `ff85983485002a5eb2986e0c61d62e94` | U |  |
| `code/r38_p3_sens.py` | `1c7697146435fcb4a4dbd5bc394c5e7a` | R | full md5 in the registration record |
| `code/r38_prep.py` | `380a80cd539a6a0d0dc395e9768d92d2` | U |  |
| `code/r38_r10a.py` | `76294c0b61fea343ef0c4ea31287264d` | R | full md5 in the registration record |
| `code/r38_s0.py` | `504a48c0eeb52f6bbb9ca43be55938ee` | R | full md5 in the registration record |
| `code/r38_s0_bias2.py` | `61229439eea8ba434e026325775aa8fd` | R | full md5 in the registration record |
| `code/r38_static.py` | `9deb1c9f96de257aeb118d13aaf247f8` | R | full md5 in the registration record |
| `code/r38_verify_b.py` | `beb68f06ad965a3c9a8abd6c4567db95` | R | full md5 in the registration record |
| `code/r38_verify_c.py` | `784aef82b0637c344534f27ef2649e39` | R | full md5 in the registration record |
| `code/r38_verify_p1.py` | `22f4e1d568e6455c82d72d98e4cc0539` | R | full md5 in the registration record |
| `code/r38_verify_p1n.py` | `1b86e96f708079c4c2916a2d3e14667c` | R | full md5 in the registration record |
| `code/r39_D_stem.py` | `2443e391c25a57af588ab6631208e784` | R | full md5 in the registration record |
| `code/r39_D_unfold.py` | `26c22d770668b0a8f905fa1551f05c2c` | R | full md5 in the registration record |
| `code/r39_core.py` | `04e57acc63848d37fb57cae3a6bb02d0` | R | full md5 in the registration record |
| `code/r39_e5prep.py` | `aee877c8bd8a5d14717cf6262b84928b` | R | full md5 in the registration record |
| `code/r39_e6.py` | `6729de931b792dd6d031ca57992e9402` | R | full md5 in the registration record |
| `code/r39_r10a_posthoc.py` | `5f20e8591412607000320be0edf7b561` | R | full md5 in the registration record |
| `code/r39_r6.py` | `a12e10a5095640c077badd83f93da713` | R | full md5 in the registration record |
| `code/r39_roundtrip.py` | `a4984f12725e509edbd908fa7ad4fff4` | R | full md5 in the registration record |
| `code/r39_sec6.py` | `1b11b675bfd25b2e34dafebc217b5618` | R | full md5 in the registration record |
| `code/r41_compute.py` | `ac8fc3ce936e13216552da71df8bcb8d` | R | full md5 in the registration record |
| `code/r41_compute_prep.py` | `add2ad62ab44ea5fc513a766a00bd7a2` | R | full md5 in the registration record |
| `code/r41_dcomp.py` | `eff667970195cb45117f89b6fc9fbca5` | R | full md5 in the registration record |
| `code/r41_dff50.py` | `88360ac232cad8a36a207a741f3bcdb9` | R | full md5 in the registration record |
| `code/r41_gapinfo.py` | `b823de6670721190025b506f3a79b475` | R | full md5 in the registration record |
| `code/r41_synth.py` | `982f6f07ac0f92d342e95b98f04cd70a` | R | full md5 in the registration record |
| `code/r41_synth_adv.py` | `96e6bc9e16cfbcc80512baa2d020fd0c` | R | full md5 in the registration record |
| `code/r43_sens.py` | `1643ce8b2bc0133a4a7cd37b4a16de60` | R | full md5 in the registration record |
| `code/r44_b6size_cell26.py` | `41fb4e7319e297133e2dcc0e16ce1a33` | U | notebook cell 26 (2026-10-28 06:04:37 UTC) verbatim below a header added for this repository |
| `code/r44_corr.py` | `4c23a80cca4644bd52a03ffc4f3b7d19` | A | addendum 04:47:45 UTC |
| `code/r44_corr2.py` | `972af7945c56590d4d8ef15f964a5ce5` | P | md5 prefix written to the run log `results/logs/r44.log` at the start of the queue, 2026-10-28 06:50:45 UTC (record addendum4 names the script) |
| `code/r44_corr_v2.py` | `f364e242d17ae36075c5c761a0d76fba` | A | addendum2 05:56:13 UTC |
| `code/r44_inv.py` | `b4562167ff541206c287f9d43282c79a` | A | addendum 04:47:45 UTC |
| `code/r44_p1comp2.py` | `06b23792e1e922755be61e70678f5c25` | U | replaces part p1comp of `r44_w2p1.py` as described in `_addendum5_fix1.json`; md5 printed at upload before the run (2026-10-28 09:02 UTC) |
| `code/r44_p1win.py` | `90b241be249fe0c2f94d1d0f8cf4a8d4` | A | addendum2 05:56:13 UTC |
| `code/r44_pairw.py` | `208d60feb7ef118952d0c245f279497c` | A | addendum2 05:56:13 UTC |
| `code/r44_redsame.py` | `2a6a20cf5687a08611518b0454174a3f` | U | md5 prefix written to `results/logs/r44.log` at the start of the queue (11:55:58 UTC); two earlier starts stopped before computing (missing package; a variable-name clash fixed in this version) |
| `code/r44_repro.py` | `a8bd5a23aa5f896d1e267dac160fedc0` | U | implementation check, run on 2026-10-28 |
| `code/r44_seqid_cell30.py` | `baba0c25848e601ce7207da39d4fa636` | U | notebook cell 30 verbatim below a header added for this repository; not registered |
| `code/r44_w2p1.py` | `70d6d7f49b3847a43f947192a4506f58` | P | md5 prefix written into the correction record `registration/prereg_2026-10-28_addendum5_fix1.json` (08:59:08 UTC, after part w2 had run at 08:55 UTC); part w2 ran as registered, part p1comp stopped before producing values |
| `code/r45_f1.py` | `55889bfa3aa916f3cabdf6d5517434b6` | A | addendum7 2026-10-28 15:39:50 UTC (`script_md5`) |
| `code/r45_folded_frac.py` | `8420f85172445aaada47480d2e7a89e8` | U | notebook run of `/content/folded_frac.py` (2026-10-29 11:31–11:37 UTC) verbatim below a header added for this repository; not registered |
| `code/r45_mpstart.py` | `642e12e591c07a3f9d907a96166f05a1` | A | addendum8 2026-10-29 11:15:02 UTC (`script_md5`) |
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
| `registration/prereg_round34.json` | `0453e32ea141d0c8d6d1258542e3ebd8` | P | checksum manifest sent with the code, 2026-10-19 |
| `registration/prereg_round37.json` | `3692c5e053c1a3fd87a22df2eefa34e7` | U |  |
| `run_records/queues/q_b6_cpu.sh` | `ab6ed598d287f688b1205816f6ec7431` | R | full md5 in the registration record |
| `run_records/queues/q_c_cpu.sh` | `a83e8fd847115e5c4338d26eed3ee592` | R | full md5 in the registration record |
| `run_records/queues/q_full_gpu10.sh` | `5dfe38f138ea1f8f221cc3e568155af4` | R | full md5 in the registration record |
| `run_records/queues/q_full_gpu12.sh` | `c4c61cf280f906821c69b7e7e36bf7af` | U |  |
| `run_records/queues/q_full_gpu13.sh` | `ccc8247223c973479a9eab26bf93fe9d` | R | full md5 in the registration record |
| `run_records/queues/q_full_gpu14.sh` | `84367c203520890f48c111d4150439a7` | R | full md5 in the registration record |
| `run_records/queues/q_full_gpu6.sh` | `d1c6079e97024d7a24eaf8ee8ba2f8e4` | R | full md5 in the registration record |
| `run_records/queues/q_full_gpu8.sh` | `3247137d47931d6a21d3dec884bec206` | R | full md5 in the registration record |
| `run_records/queues/q_full_gpu9.sh` | `f2cd9db4e000eaa36bef403178ca751f` | R | full md5 in the registration record |
| `run_records/queues/q_p1_cpu.sh` | `98dafc7db4669fa0c867363d569e51a6` | R | full md5 in the registration record |
| `run_records/queues/q_p3s_cpu.sh` | `dc9952e7a634a33e7d9489cbb5644e0c` | R | full md5 in the registration record |
| `run_records/queues/q_r39_cpu1.sh` | `138ad916e3faa270352f6062a6559f40` | R | full md5 in the registration record |
| `run_records/queues/q_r39_cpu2.sh` | `a5993779ffc78a505723288a2770853c` | R | full md5 in the registration record |
| `run_records/queues/q_r39_gpu.sh` | `9a0e52c775f86001fab404835bf97dca` | R | full md5 in the registration record |
| `run_records/queues/q_r44_cpu.sh` | `7c43fc43d4263cdfb24a65da8b02419e` | U |  |
| `run_records/queues/q_r44_gpu.sh` | `c535049adaa9335d4e4422be543dc394` | U |  |
| `run_records/queues/q_r44b_corr.sh` | `12e058930ddb2843abdc0b856fcb9ece` | U |  |
| `run_records/queues/q_r44b_misc.sh` | `7112c17b8c2e46a5a17b29b8440f5cbe` | U |  |
| `run_records/queues/q_r44c.sh` | `0a582e7eab7640bb6b9e67fad46543df` | U |  |
| `run_records/queues/q_r44d.sh` | `c36712c2119c459dcfaa25a13eb0ced5` | U |  |
| `run_records/queues/q_r45.sh` | `3b1de9082ceff5a5b51ce9f6cd2071da` | A | addendum7 (`queue_md5`) |
| `run_records/queues/r39_resume_20261021.sh` | `a32e0766f12b37572dbc215c9dbe295c` | U |  |
