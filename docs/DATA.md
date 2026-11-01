# Data

Three kinds of data are involved: public inputs that anyone can download from the original source,
intermediate tables that were made during the project and are needed to rerun the scripts, and the outputs
behind the numbers in the manuscript. Paths are relative to the project root `ROOT` (see the README).

Status column: **public** = download from the source; **deposit** = included in the data deposit (Zenodo,
DOI to be assigned at submission); **regenerate** = written by a script in this repository from the items
above.

## 1. Public inputs

| Data | Used for | Path under ROOT | Source and license |
|---|---|---|---|
| mdCATH (CHARMM22*/TIP3P, 320/348/379/413/450 K, 5 replicas per domain) | all same-sequence temperature cells, homologue pairs, P1, P1′, C | `raw/mdcath/mdcath_dataset_<domain>.h5`; listing `raw/mdcath_all_domains_listing.json` | Hugging Face dataset `compsciencelab/mdCATH`, CC BY 4.0. `r38_c_fetch.py` and `r38_p1n_fetch.py` read only the needed frames over HTTP into `raw/mdcath_c/`, `raw/mdcath_p1n/` |
| ATLAS (CHARMM36m/TIP3P, 300 K, 3 replicas) | R3′ force-field cell, E5, AlphaFlow test set (82 chains) for the 7b box and the model-error ladder row | `raw/atlas/<chain>_protein.zip`, `raw/atlas_af82/<chain>_analysis.zip` | https://www.dsimb.inserm.fr/ATLAS, CC BY-NC 4.0 |
| AlphaFlow / ESMFlow public samples and weights | model-error ladder row, 7b box, P2 secondary comparisons | `raw/models/alphaflow/samples/<model>.zip`, `raw/models/alphaflow/params/alphaflow_md_distilled_202402.pt` | AlphaFlow repository (commit `0408d7c`) and its released samples and weights, MIT |
| aSAMt weights and mdCATH input bundle | P2 (aSAMt ensembles at 320–450 K) | `raw/models/asam/input_mdcath.zip`; code `sam2` commit `ce076ed` | aSAM (Janson et al.); repository `giacomo-janson/sam2` |
| BioEmu CATH MD (amber99sb-ildn, 300 K; `MSR_cath2`, `ONE_cath1`) | force-field + protocol ladder row, R3″ protocol-only cell | `raw/bioemu_cath/MSR_cath2.zip`, `raw/bioemu_cath/ONE_cath1.zip` | BioEmu MD data release, CATH part: Zenodo record 15629740 (doi:10.5281/zenodo.15629740), CDLA-Permissive-2.0 |
| MegaSim (ff99SB-disp, 295 K; wild types and substitution variants) | substitution ladder row (C4), R5′ | `raw/megasim/MSR_megasim_merge.zip`, `raw/megasim/MSR_megasim_mutants_disp_allatom.zip` | BioEmu MD data release, MegaSim part: Zenodo record 15641184 (doi:10.5281/zenodo.15641184), CDLA-Permissive-2.0 |
| Tsuboyama et al. 2023 mega-scale stability (ΔG, ΔΔG) | C4 (ΔΔG of MegaSim variants) | `raw/tsuboyama/Processed_K50_dG_datasets.zip` | Zenodo record 7992926, CC BY 4.0 |
| MFM all-atom reference MD (ff99SB-disp/TIP3P, 150 mM, 300 K; "SM" in the ladder) | force-field + protocol ladder row | `raw/mfm_sm/<domain>.npz` (positions, atom lines, metadata) | transferable-cg data release (Park et al.; Hugging Face `abpark/transferable-cg`), MIT. The conversion of the released files into these npz files was done in a notebook cell; the 54-domain npz set is deposited |
| Majewski fast-folder MD (CHARMM22*/TIP3P, 350 K) | R7 unfolding (upper anchor) | `raw/majewski/` | torchmd-protein-thermodynamics data release (see its `Datasets/README.md`), MIT |
| CATH domain lists and boundaries | superfamily labels, chain-linking candidates | `raw/cath/cath-domain-list-v4_2_0.txt`, `raw/cath/cath-domain-list.txt`, `raw/chaining/cath-domain-boundaries.txt` | CATH (v4.2.0 list as used by mdCATH; the unversioned list and the CathDomall v4.4 boundaries file as downloaded in 2026-10). The deposit manifest records their md5 so the same files can be identified |
| RCSB PDB entries (150 two-domain chains) | chain linking (supplement) | `raw/chaining/pdb/<id>.pdb.gz` | downloaded by `chaining.py` from files.rcsb.org |
| SIFTS `pdb_chain_taxonomy.tsv.gz` | organism mapping for P1′ pairs | `raw/` (4 MB) | EBI SIFTS |
| TEMPURA; Engqvist 2018 growth temperatures | optimal growth temperature of source organisms | used to build `results/pilot_mdcath_organism_temp_v2.tsv` | TEMPURA database; Engqvist 2018 (Zenodo record 1175609) |

## 2. Intermediate tables to deposit

These were made in notebook cells or by earlier code that is not part of this repository, and the scripts
read them as inputs. Deposit them as they are; their checksums go into the deposit manifest.

| File | What it is | Read by |
|---|---|---|
| `results/mdcath_5T_setup.json` | mdCATH domain list (336 domains), temperatures, frame step | `mdcath_5T_*`, `c1_core.py` and everything built on `derived/mdcath_5T` |
| `results/mdcath_spread_priority.json` | processing order for `mdcath_spread*.py` | `mdcath_spread.py`, `mdcath_spread_hot.py` |
| `results/pilot_mdcath_homolog_pairs.tsv` | same-superfamily mdCATH domain pairs (H-pool) | `p2_gap_consensus.py`, `r38_p3.py`, `r38_c_probe.py` |
| `results/pilot_mdcath_organism_temp_v2.tsv` | source organism, taxid, growth temperature per domain | `r38_p1.py`, `ht_theta_meso.py` |
| `results/overlap_mdcath_atlas.json` | mdCATH ∩ ATLAS chains with sequence coverage | `r37_ff.py`, `r38_dff*.py`, `r39_e5prep.py`, `r41_dff50.py` |
| `results/ht_shared_spread.tsv`, `results/ht_theta_star.tsv` | H-T thermophile–mesophile pairs and first θ* estimates (2026-10-17) | `ht_theta_meso.py`, `folded_analyses.py`, `ht_own_alignment.py` |
| `results/r37_trunk_candidates.tsv`, `results/r37_trunk_pairs.json`, `results/r37_trunk_pairs.tsv` | trunk pair selection (round37) | `r38_static.py`, `r38_b6.py`, `r38_r10a.py`, `r39_sec6.py` |
| `raw/asamt_test90_needs.json` | aSAMt test-set domains and temperatures needed | `r38_prep.py`, `r38_md_ref_p90.py` |
| `derived/r38/*.txt` (`p2_67`, `asamt_sanity5`, `q3_p2_67_ordered`, `q5a_450_17`, `q5b_p90_4`, `rest18`, `q_cond_450_50`, `p2_sf23`, `p2_no_big4`, `p2_p80`), `derived/r38/alphaflow_input.csv` | domain lists for the P2 generation and aggregation queues (the 67 P2 domains, the 5 correctness-check domains drawn with seed 20261019, subsets for sensitivity versions) | `r38_s0.py`, `r38_asamt_gen.py`, `r38_alphaflow_gen.py`, `r38_p2_*`, `r39_core.py`, `r39_e6.py`, `r43_sens.py` |
| `derived/r38/p1n/pairs.json` | P1′ pairs, recorded 2026-10-20 13:35:54 UTC | `r38_p1n_fetch.py`, `r38_p1n.py` |
| `results/v_residue_unit_graded.tsv`, `results/xprot_bench.tsv`, `derived/xprot_structs.npz`, `derived/atlas_pdb/` | ruler-benchmark pairs and structures (1,674 pairs; 208-structure block) | `hu1_rulers.py`, `dfrechet_1674.py`, `analysis_fold208*.py` |
| `results/pu_fold208_setup.json`, `derived/pu_fold208_ca.npz`, `results/pu_tm_factorial_setup.json`, `results/pu_tm_factorial_B_setup.json`, `derived/pu_blockB_ca.npz` | structure blocks for the TM factorial design and fold retrieval | `pu_*.py` |
| `raw/mfm_sm/*.npz` | see section 1 | `ladder_ff*.py` |
| `results/ladder_*_summary.json` | aggregated ladder rows (aggregation cells not kept) | text of result 4.8 |
| `results/prereg_2026-10-17.json` | the first time-stamped record (all rounds, amendments, checklists) | `make_s1.py` (through `entries_timeline.tsv`) |
| `results/prereg_2026-10-28_addendum.json`, `_addendum2.json` … `_addendum6.json`, `_addendum5_fix1.json`; `results/prereg_2026-10-29_addendum7.json`, `_addendum8.json` | the six round44 records and one correction (2026-10-28) and the two round45 records (written 2026-10-28 15:39:50 and 2026-10-29 11:15:02 UTC), separate files written after the first record was closed; copies in `registration/` | `make_s1.py` |
| `results/prereg_2026-10-30_addendum9.json`, `_addendum9_fix1.json`, `_addendum10.json`, `_addendum10_fix1.json`, `_addendum10_fix2.json`, `_addendum11.json` | the three round46–48 records and their three corrections (written 2026-10-29 16:34:05 to 2026-10-30 07:30:01 UTC); the file-name date of the round45–48 records is KST; copies in `registration/` | `make_s1.py` |
| `results/prereg_2026-10-31_addendum12.json`, `_addendum12_fix1.json` | the round49 record and its correction (written 2026-10-31 13:51:43 and 14:25:18 UTC); copies in `registration/` | `make_s1.py` |
| `results/_registration_evidence_2026-10-26/` | evidence for the times of the first record: its 101 retained Google Drive versions (zip) with their manifest, the manifest of its 46 snapshot copies, and the table of the first appearance and later changes of every entry (`entries_timeline.tsv`) | `make_s1.py` |

## 3. Outputs

| What | Path | Status |
|---|---|---|
| Summary tables behind every number in the text | `results/*.json`, `results/*.tsv`, `results/*.txt` (map in `docs/RESULTS_MAP.md`) | deposit (small) |
| Backups written next to each result at run time | `results/_backup/<date>/` | deposit |
| Run logs and completion stamps | `results/logs/*.log`, `results/logs/*.done` | deposit |
| Generated model ensembles (aSAMt 320–450 K, AlphaFlow distilled) | `derived/r38/asamt/`, `derived/r38/alphaflow/` | deposit (stochastic; regenerating them gives different samples) |
| Per-domain statistics (JSON) | `derived/r38/*/`, `derived/r39/*/`, `derived/r41/*/`, `derived/r44/*/`, `derived/r45/*/`, `derived/r46/`, `derived/r47/`, `derived/r48/`, `derived/r49/`, `derived/r37_*/`, `derived/ladder_*/`, `derived/c4_sub/`, `derived/c1_core/` | deposit the JSON files; regenerate the rest |
| Extracted mdCATH frames and spreads (npz) | `derived/mdcath_5T/`, `derived/mdcath_5T_folded/`, `derived/mdcath_spread*/` | regenerate from mdCATH (large) |

`tools/pack_for_zenodo.py` builds the deposit from a project root following this table and writes a checksum
manifest. Run it on the machine that holds the project directory (Colab with Drive mounted, in our case).
