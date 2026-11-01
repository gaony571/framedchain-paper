# Environments

All runs used Google Colab with the project root on Google Drive (`/content/drive/MyDrive/atlas_framed_curves`).
GPU stages ran on A100-SXM4-40GB runtimes (12 vCPUs, Intel Xeon @ 2.20 GHz); the session outputs give the GPU
model on 19, 20, 21, 22 and 28 October 2026 and the driver, 580.82.07, on 21 and 22 October. CPU-only scripts ran
on whichever runtime was attached, including CPU-only runtimes with 2 vCPUs (for example on 16 October, with
torch 2.11.0+cpu) and the CPUs of an L4 runtime (round43, 22 October; 12 vCPUs, Intel Xeon @ 2.20 GHz). The
equipment check of the compute table is in `docs/INTERACTIVE_CHECKS.md` 3.

## 1. Base Colab environment (most scripts)

Packages used (with the script that needs each optional one): `requirements-base.txt`. In addition the scripts
import the in-repository `framespace_pkg` (framespace 0.2.0) and `framespace_wip` through `sys.path`.
`requirements-check.txt` lists the versions with which the tests and `examples/demo_synthetic.py` were checked
on 2026-10-26 (CPU only).

Versions recorded during the runs:

| Package | Version | Where recorded |
|---|---|---|
| numpy | 2.1 | 7b box environment check, 2026-10-19 |
| numpy | 2.1.3 | equipment check of the compute table, 2026-10-22 (`docs/INTERACTIVE_CHECKS.md` 3) |
| mdtraj | 1.11.2 | 7b box environment check, 2026-10-19 |
| torch | 2.11.0+cu130 | compute table (round41 g5), 2026-10-22 |
| h5py | 3.14.0 | equipment check of the compute table, 2026-10-22 (`docs/INTERACTIVE_CHECKS.md` 3) |
| tmtools | 0.3.0 | project record |

A full `pip freeze` of the base environment was **not** recorded at run time. The resume script
(`run_records/queues/r39_resume_20261021.sh`) printed a `BASECHECK` line with the mdtraj, numpy and torch
versions to `results/logs/r39_resume.log`, which is part of the data deposit. For a fresh installation,
`tools/record_env.py` writes the same kind of record for whatever environment you use.

The 7b box showed why this matters: the AlphaFlow analysis code gives different model-side metrics under
numpy 2.1 / mdtraj 1.11.2 than under its own README environment (mdtraj 1.9.9), because `atom_slice` keeps a
different atom order (results section 5.1). The standard-metric runs were therefore done in afenv.

## 2. afenv — AlphaFlow generation and the AlphaFlow analysis code

- Python 3.9, torch 1.12.1+cu113, numpy 1.21.2, scipy 1.7.1, mdtraj 1.9.9, pandas 1.5.3, pytorch-lightning
  2.0.4, torchmetrics 1.0.3, OpenFold at commit `103d037` built from source with its CUDA kernels.
- `afenv_freeze.txt` — `pip freeze` of the environment as used.
- `setup_afenv3.sh` — the third and final build step (OpenFold rebuilt with gcc-11 because nvcc 11.8 does not
  accept the runtime's gcc 13). The earlier steps that created the environment under `/content/afenv` were run
  from notebook cells; `afenv_freeze.txt` is the authoritative package list.
- `af_bar_readmeenv.sh` — reruns the AlphaFlow table (`scripts.analyze_ensembles`, `scripts.print_analysis` at
  AlphaFlow commit `0408d7c`) in this environment (results section 5.1).
- Used by `r38_alphaflow_gen.py` and `r38_7b_split.py` (they call `/content/afenv/bin/python`) and by
  `af_bar_readmeenv.sh`. `r38_af_bar.py` itself ran in the base environment; comparing its output with the
  afenv rerun is how the environment effect in section 5.1 was found.

## 3. samenv — aSAMt generation

- Python 3.12 (the sam2 code at `ce076ed` needs PEP 701 f-strings), torch 2.6.0+cu124, numpy 2.5.4,
  mdtraj 1.11.2, diffusers 0.41.0, sam2 installed editable from commit `ce076ed`.
- `samenv_freeze.txt` — `pip freeze` as used. `setup_samenv3.sh` — build script (micromamba, conda-forge).
- Used by `r38_asamt_gen.py` (calls `/content/samenv/bin/python` on `generate_ensemble.py`).

## 4. Other tools

| Tool | Version | Used for |
|---|---|---|
| AlphaFlow | commit `0408d7c`, cloned to `ROOT/code/ext/alphaflow` | generation, analysis code of the 7b box |
| sam2 (aSAM / aSAMt) | commit `ce076ed` | aSAMt ensembles |
| OpenFold | commit `103d037` | AlphaFlow dependency |
| US-align | commit `23f14ff` (built 2026-10-19), `/content/usalign/USalign` | timing comparison in `r39_sec6.py` and `r41_compute.py` |
| Foldseek | `10-941cd33` | timing comparison in `r41_compute.py` |
| tmtools | 0.3.0 | TM-align alignments inside Python scripts |
| Node.js | any recent version (checked with v24.19.0) | `code/negtype/*.js`, run from that directory: `node so3fourier.js`, `node len2check.js`, `node negtype.js [search maxLen N trials]` (seconds each) |

No container image is provided: the two model environments need CUDA builds (OpenFold kernels, torch with
cu113/cu124) and model weights whose redistribution is governed by their own terms. The setup scripts, freeze
files and `requirements-*.txt` document what was installed.
