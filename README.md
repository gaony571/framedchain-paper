# framedchain-paper

Analysis code, run records and registration records for the paper

**Comparing protein conformational ensembles across sequences and lengths**
(Hyeonsu An; in preparation)

which represents protein backbones of any length as framed chains, builds a metric space and ensemble statistics on
them, and uses them to compare structural ensembles across temperatures, force fields, generative models and
sequences. The abstract is at the end of this file.

Version 1.0.0 — the code as run for the manuscript. Author: Hyeonsu An. Licence: MIT.

This repository is the **as-run** code: every script here is the version that produced the numbers in the
manuscript, byte for byte where a checksum was registered (see `docs/PROVENANCE.md`). When a registered script
had to be corrected, the corrected version and its checksum were registered as a separate, time-stamped entry;
analyses added after a result are separate files and say so in their docstrings. The reusable, documented and
tested implementation of the method is the separate package
[framedchain](https://github.com/gaony571/framedchain) (version 0.4.1); start there if you want to apply the
method to your own data.

**Names.** The package was called framespace up to version 0.3.3 and this repository framespace-paper; both were
renamed on 2026-10-26 (GitHub redirects the old addresses). The as-run copies keep the old name
(`framespace_pkg/framespace/` is framespace 0.2.0, the scripts `import framespace`, and the research module is
`framespace_wip/`), because their bytes are registered.

**Language of the comments.** Docstrings and comments inside the scripts are in Korean, the working language
of the project; they are not translated in place because the files are kept byte-for-byte as run.
[`docs/SCRIPTS.md`](docs/SCRIPTS.md) describes every script in English (purpose, stages, inputs, outputs,
registration entry).

## Quick start (no downloads, CPU, a few minutes)

```bash
pip install -r environment/requirements-check.txt
python examples/demo_synthetic.py            # compare with examples/demo_expected_output.txt
PYTHONPATH=framespace_pkg python -m pytest -q framespace_wip/tests_stage0
```

The pinned versions need Python 3.12 or newer. On Windows set the path first (`set PYTHONPATH=framespace_pkg` in
cmd, `$env:PYTHONPATH = "framespace_pkg"` in PowerShell) and run `python -m pytest ...` without the prefix.

- `examples/demo_synthetic.py` applies the paper's computations to small synthetic ensembles using only the
  as-run libraries: invariances and triangle inequality of D, the reduced computation (identical to
  `red_pairs` in `code/r38_s0.py`) against the full search, and the energy-distance statistic Δ² for a null
  comparison, a shift in fluctuations and an insertion of four edges. It is deterministic and takes about
  30 s; its output on the checked versions is in `examples/demo_expected_output.txt`.
- The 25 tests of the research module run against the as-run framespace 0.2.0 (about 3 minutes).
- Both were checked on 2026-10-26 with Python 3.12.13 and the versions in
  `environment/requirements-check.txt`. The method-level tests (definitions, the reduced computation, the
  ensemble statistics fixed for this paper) are in the framedchain package (`tests/`, 46 tests in 0.4.1).

## What is here

| Path | Content |
|---|---|
| `code/` | Analysis scripts (Python), one per experiment or stage; `code/negtype/` holds three small Node.js checks. English descriptions in `docs/SCRIPTS.md` |
| `framespace_pkg/framespace/` | framespace 0.2.0 exactly as used in the runs (tag `v0.2.0`, commit `866995c` of the framedchain repository; bundle md5 `4edba596…`) |
| `framespace_wip/stage0/`, `framespace_wip/tests_stage0/` | Research module used by the runs (reduced computation, rulers, ESA/W2 helpers) and its tests (bundle md5 `6703fa26…`) |
| `framespace_wip/experiments/` | Constructed counterexample for the Sippl linear form |
| `examples/` | Synthetic demonstration and its expected output |
| `environment/` | Package lists (`requirements-base.txt` for the analysis scripts, `requirements-check.txt` for the demo and tests), setup scripts and `pip freeze` records of the two model environments; see `environment/README.md` |
| `registration/` | Local copies of registration entries (round34, round37). The full record is `results/prereg_2026-10-17.json` in the data deposit |
| `run_records/queues/` | The shell queues that ran the scripts on Colab, with the exact arguments |
| `docs/RESULTS_MAP.md` | Manuscript section → scripts → output files → registration entry |
| `docs/SCRIPTS.md` | What every script computes, in English |
| `docs/DATA.md` | Input data, intermediate tables and outputs; what is deposited and where |
| `docs/PROVENANCE.md` | Per-file checksum and registration status |
| `docs/STATEMENTS.md` | Draft Data Availability and Code Availability statements and the software references |
| `tools/` | `pack_for_zenodo.py` (bundle results and selected intermediates for the data deposit), `record_env.py` (record an environment) |
| `MANIFEST.md5` | md5 of every file in this repository (`md5sum -c MANIFEST.md5`) |

## How the code is laid out and run

All scripts take the project root as their first argument and read and write under it:

```
ROOT/
  code/              this repository's code/ (scripts load each other from here with exec)
  code/ext/alphaflow AlphaFlow repository, commit 0408d7c (only for the AlphaFlow steps)
  framespace_pkg/    this repository's framespace_pkg/
  framespace_wip/    this repository's framespace_wip/
  raw/               downloaded inputs (docs/DATA.md)
  derived/           per-domain intermediates written by the scripts
  results/           summary tables, registration record, logs (results/logs/*.done)
```

A typical call is

```bash
python code/r38_p3.py "$ROOT" pairs,dom,pibar,pair,stats 3
```

where the second argument lists stages and the optional third is the number of worker processes. Stages are
resumable: finished per-domain files are skipped. Scripts add `ROOT/framespace_pkg` and `ROOT/framespace_wip`
to `sys.path` themselves, so do not install another framespace version (framedchain 0.3.x or earlier) into the same environment when reproducing the
runs. Several scripts load the functions of another script with `exec` (for example `r38_b6.py` loads the first
part of `r38_s0.py`, and `r38_p1.py` loads `r38_p3.py`); keep the `code/` directory flat and unchanged.

### Order of the main stages

1. **Inputs.** Download the public data into `raw/` (docs/DATA.md) and place the deposited intermediate tables
   in `results/` and `derived/`.
2. **mdCATH extraction** (CPU): `mdcath_5T_extract.py` → `mdcath_5T_folded.py` → `c1_core.py`;
   `mdcath_spread.py` and `mdcath_spread_hot.py` for the per-edge spreads.
3. **Ladder and early checks** (CPU): `ladder_*`, `c4_substitution.py`, `e0_resolution.py`, `ladder_half.py`,
   `ht_theta_meso.py`, the ruler comparisons (`hu1_*`, `pu_*`, `dfrechet_1674.py`, `chaining.py`).
4. **Round37 CPU ladder:** `r37_trunk.py`, `r37_ff.py`, `r37_protocol.py`, `r37_subst.py`, `r37_unfold.py`.
5. **Round36 GPU plan:** first session `r38_s0.py`; model ensembles `r38_prep.py`, `r38_md_ref*.py`,
   `r38_asamt_gen.py`, `r38_alphaflow_gen.py`; comparisons `r38_compare*.py`, `r38_p2_*`; primary experiments
   `r38_p3.py`, `r38_p1.py`, `r38_c.py`, `r38_p1n.py`, `r38_b6.py`, `r38_static.py`, `r38_r10a.py`; the
   standard-metric box `r38_af_bar.py`, `r38_7b_*`. P3/B6, P1, C and P1′ each have a `r38_verify_*.py` that
   checks the implementation on synthetic data with known answers before the registered run.
6. **Round39 scale experiments:** `r39_e5prep.py`, `r39_core.py`, `r39_D_stem.py`, `r39_D_unfold.py`,
   `r39_e6.py`, `r39_r6.py`, `r39_sec6.py`, `r39_roundtrip.py`.
7. **Round41 and round43:** `r41_*`, `r43_sens.py`.

The exact sequences that ran are in `run_records/queues/` (`$R` = ROOT, `$C` = ROOT/code, `$L` = log directory).

## Environments and hardware

The runs used Google Colab (A100-SXM4-40GB for the GPU stages, 12 vCPU Intel Xeon) with the project directory
on Google Drive. Three Python environments were used:

- **Base Colab environment** (Python 3.12) for almost all scripts; packages in
  `environment/requirements-base.txt`. Recorded versions: numpy 2.1 and mdtraj 1.11.2 (2026-10-19),
  torch 2.11.0+cu130 (2026-10-22), tmtools 0.3.0. A full `pip freeze` of this environment was not recorded at
  run time; `tools/record_env.py` writes one for the environment you use.
- **afenv** (Python 3.9, torch 1.12.1+cu113, OpenFold `103d037`) for AlphaFlow generation and the AlphaFlow
  analysis code: `environment/afenv_freeze.txt`, `environment/setup_afenv3.sh`.
- **samenv** (Python 3.12, torch 2.6.0+cu124, sam2 `ce076ed`) for aSAMt generation:
  `environment/samenv_freeze.txt`, `environment/setup_samenv3.sh`.

External tools: AlphaFlow commit `0408d7c`, sam2 (aSAM) commit `ce076ed`, OpenFold `103d037`, US-align
`23f14ff`, Foldseek `10-941cd33`, tmtools 0.3.0. See `environment/README.md`.

Rough cost on that hardware: the CPU stages take minutes to hours each; the GPU stages (model generation,
reduced computation over millions of frame pairs) took days in total. The compute table of the paper
(`results/r41_compute.json`, result 2.6) gives per-pair timings.

## Reproducing at three levels

1. **Check the numbers in the text** from the deposited `results/` tables: each row of `docs/RESULTS_MAP.md`
   names the file that holds the value.
2. **Re-run a statistic from deposited intermediates** (CPU, minutes to hours): `stats` stages read the
   per-domain files in `derived/`, for example `python code/r38_p3.py "$ROOT" stats`.
3. **Re-run from raw data** (GPU, days): download the inputs and run the stages in order. The model ensembles
   (aSAMt, AlphaFlow) are stochastic; the deposit contains the generated ensembles that were analysed.

## Known gaps

- Some inputs were produced in notebook cells that were not kept as scripts. They are deposited as data and
  listed in `docs/DATA.md` (domain and pair lists, the mdCATH setup file, the ruler-benchmark structure sets,
  the growth-temperature table, the summary files of the old ladder rows).
- Two early checks (TM-align comparison on 12 pairs, reversed chains) were computed interactively with
  `framespace_wip/stage0/rulers.py`; their tables are in the project record only.
- The base Colab environment was not frozen at run time (see above).
- Some paths are Colab-specific (`/content/...` scratch directories, `/usr/bin/python3`, the Drive project path
  in `environment/af_bar_readmeenv.sh`). They are working directories, not inputs; change them if you run
  elsewhere.
- Scripts whose checksum was not registered are marked in `docs/PROVENANCE.md`.

## Licence and citation

MIT licence, copyright (c) 2026 Hyeonsu An (`LICENSE`). Third-party code is not included: AlphaFlow, sam2,
OpenFold, US-align and Foldseek are obtained from their own repositories under their own licences. If you use
this code, please cite the accompanying paper (in preparation) and the software (`CITATION.cff`; Zenodo DOI to
be added at release).

## Abstract

Proteins are increasingly described as conformational ensembles, but no common geometry holds structures and ensembles of different sequences and lengths: most ensemble measures assume identical sequences, and residue-pair comparisons discard residues without a partner. Here I represent protein backbones of any length as framed chains. Each virtual Cα–Cα bond becomes a 3×3 matrix, the orthonormal frame of the bond and its peptide plane scaled by bond length; an absent bond is the zero matrix, the only rotation-invariant choice. Minimizing Euclidean distance between chains over global rotations and insertions of absent bonds gives an order-preserving edit-distance metric across lengths; at any optimal alignment, which can be certified, its square splits exactly into shared-bond and unshared-bond parts. Likewise, the score that TM-align maximizes, with its distance scale fixed, converts exactly into an edit-distance metric on α-carbon positions modulo rigid motion. For ensembles, energy statistics corrected for replicate noise are unbiased under any underlying metric. Calibrated on same-sequence temperature responses, they prove invariant to replica and frame numbers but grow with domain size, which comparisons must account for. Once correspondence and rotation are fixed, within or across sequences, the space is Euclidean, so mean shifts and fluctuations can be compared bond by bond. On held-out domains, the mean-shift errors of a temperature-conditioned generative model match or exceed the same domains’ shifts under 59 K of heating, yet point largely away from the heating direction. In a pre-specified replication across 32 superfamilies, thermophilic domains fluctuate less on shared bonds than mesophilic homologs at equal temperature. In a confirmatory set of 4,152 pairs of larger homologous domains, shared-bond fluctuation profiles grow modestly more similar with sequence identity beyond shared-core structure, elastic-network predictions and nonlinear structural covariates; in shorter domains, structure explains the association. Ensembles of different sequences thus share one geometry.
