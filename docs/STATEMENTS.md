# Availability statements and software references (PLOS Computational Biology)

## What the journal requires (code availability policy, checked 2026-10-26)

Source: https://journals.plos.org/ploscompbiol/s/code-availability

| Requirement | How it is met |
|---|---|
| Code that directly relates to the results is shared upon publication; editors and reviewers can access it at submission (a private link is acceptable) | This repository (all analysis scripts as run) and the framedchain package; reviewer access at submission, public at publication |
| The Data Availability Statement says where and how the code can be accessed, with a repository name and a DOI or other persistent identifier; placeholder text is allowed if the DOI is issued at publication | Draft below; Zenodo DOIs to be reserved before submission (placeholders marked XXXXXXX) |
| A licence is clearly specified; an OSI-conforming open-source licence is strongly encouraged; Creative Commons is not recommended for code | MIT (`LICENSE`) for both repositories |
| A repository that issues DOIs and supports versioning is recommended (GitHub alone has no DOI; pair it with Zenodo); cite the code in the reference list | GitHub releases archived on Zenodo; references below |
| Dependencies, run environment, versions and system requirements are documented; a standard test dataset is included where feasible; code is clearly annotated | `README.md`, `environment/`, `examples/demo_synthetic.py` (synthetic, no downloads), `docs/SCRIPTS.md` (English description of every script; in-file comments are Korean) |
| Preprocessing code is strongly encouraged | Data-preparation scripts are included; inputs made in notebook cells are deposited as data (`docs/DATA.md`) |

## Data Availability Statement (draft)

> All input data are publicly available from their original sources: mdCATH (Hugging Face dataset
> compsciencelab/mdCATH), ATLAS (https://www.dsimb.inserm.fr/ATLAS), the BioEmu molecular dynamics data
> (Zenodo records 15629740 and 15641184), the Tsuboyama et al. mega-scale stability data
> (Zenodo record 7992926), the transferable-cg reference simulations (Hugging Face dataset
> abpark/transferable-cg), the torchmd-protein-thermodynamics fast-folder simulations, the AlphaFlow and aSAM
> released samples and weights, CATH, SIFTS and the RCSB PDB. All data generated in this study that are
> needed to reproduce the reported values — summary tables for every result, per-domain statistics, the
> generated aSAMt and AlphaFlow ensembles, the domain and pair lists, and the time-stamped pre-registration
> record — are deposited at Zenodo (DOI: [10.5281/zenodo.XXXXXXX]). The analysis code used to produce all
> results is available at https://github.com/gaony571/framedchain-paper and archived at Zenodo (DOI:
> [10.5281/zenodo.XXXXXXX]); the method is implemented in the framedchain Python package
> (https://github.com/gaony571/framedchain, version 0.4.1, archived at Zenodo, DOI: [10.5281/zenodo.XXXXXXX]).
> Both are released under the MIT licence. A file-level map from each result to its script and data file is
> given in docs/RESULTS_MAP.md of the code repository.

## Code availability sentence for the Methods (draft)

> The analysis scripts exactly as run, with the run queues, environment records and checksum records, are
> available at https://github.com/gaony571/framedchain-paper (version 1.2.0; MIT licence; Zenodo DOI:
> [10.5281/zenodo.XXXXXXX]). The reusable implementation of the space, the distance and the ensemble
> statistics is the framedchain package (version 0.4.1; MIT licence; Zenodo DOI: [10.5281/zenodo.XXXXXXX]).

## Software references (for the reference list)

> An H. framedchain: framed chains for protein structures and ensembles of any length [software]. Version 0.4.1.
> Zenodo; 2026. doi:[10.5281/zenodo.XXXXXXX]

> An H. framedchain-paper: analysis code for "Comparing protein conformational ensembles across sequences and lengths" [software].
> Version 1.2.0.
> Zenodo; 2026. doi:[10.5281/zenodo.XXXXXXX]

> An H. Data for "Comparing protein conformational ensembles across sequences and lengths" [dataset].
> Zenodo; 2026. doi:[10.5281/zenodo.XXXXXXX]

## Status and remaining steps

Done (2026-10-26): both repositories are on GitHub as private repositories with version tags —
`gaony571/framedchain` (tags `v0.2.0` = the version used in the runs, `v0.3.0` to `v0.4.1`) and
`gaony571/framedchain-paper` (tags `v1.0.0` = first push, `v1.1.0` = with the round44 checks of 2026-10-28, `v1.2.0` = with the round45 checks of 2026-10-28 and 10-29). Both were renamed from framespace and framespace-paper on
2026-10-26 (GitHub redirects the old addresses). No GitHub *release* has been published yet, on purpose: Zenodo
archives a repository only for releases published after the repository is connected.

Remaining, after the manuscript is written and posted (these need the paper or the DOIs):

1. Connect both repositories to Zenodo, then publish a GitHub release of each from its latest `main`
   (framedchain: tag `v0.4.1`, or a later one if the package changes; framedchain-paper: `v1.2.0`, or, if `main` changes after 2026-10-29, a tag of the then-current `main`,
   e.g. `v1.2.1`, with the version updated in `CITATION.cff` and the README).
   Zenodo mints a DOI for each release.
2. Run `tools/pack_for_zenodo.py` on the project root (Colab with Drive mounted; first with `--dry-run` to see
   the size against the Zenodo limit), check `MANIFEST_data.tsv`, and upload the data deposit. ATLAS files
   (CC BY-NC 4.0) are not redistributed; the deposit lists them by name and checksum only.
3. Fill the DOIs in the statements above (the paper title is filled in), in both `CITATION.cff` files and in the READMEs.
4. At submission give the editors and reviewers access (a private link is acceptable), or make the
   repositories public when the preprint is posted.
5. Decide whether the registration record is published with the data deposit or registered separately (for
   example on OSF), and cite it in the Methods.
