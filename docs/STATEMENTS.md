# Availability statements and software references (PLOS Computational Biology)

## What the journal requires (code availability policy, checked 2026-10-26)

Source: https://journals.plos.org/ploscompbiol/s/code-availability

| Requirement | How it is met |
|---|---|
| Code that directly relates to the results is shared upon publication; editors and reviewers can access it at submission (a private link is acceptable) | This repository (all analysis scripts as run) and the framedchain package; reviewer access at submission, public at publication |
| The Data Availability Statement says where and how the code can be accessed, with a repository name and a DOI or other persistent identifier; placeholder text is allowed if the DOI is issued at publication | Draft below; the data DOI can be reserved before submission (Zenodo draft), the two software DOIs are minted when the repositories are public and released (placeholders marked XXXXXXX) |
| A license is clearly specified; an OSI-conforming open-source license is strongly encouraged; Creative Commons is not recommended for code | MIT (`LICENSE`) for both repositories |
| A repository that issues DOIs and supports versioning is recommended (GitHub alone has no DOI; pair it with Zenodo); cite the code in the reference list | GitHub releases archived on Zenodo; references below |
| Dependencies, run environment, versions and system requirements are documented; a standard test dataset is included where feasible; code is clearly annotated | `README.md`, `environment/`, `examples/demo_synthetic.py` (synthetic, no downloads), `docs/SCRIPTS.md` (English description of every script; in-file comments are Korean) |
| Preprocessing code is strongly encouraged | Data-preparation scripts are included; inputs made in notebook cells are deposited as data (`docs/DATA.md`) |

## Data Availability Statement (draft)

> All input data are publicly available from their original sources: the mdCATH simulations (Hugging Face
> dataset compsciencelab/mdCATH), the ATLAS simulations (https://www.dsimb.inserm.fr/ATLAS), CATH domain
> lists and boundaries (https://www.cathdb.info), RCSB PDB entries for the two-domain chains of the chaining test
> (https://www.rcsb.org), SIFTS organism mappings (https://www.ebi.ac.uk/pdbe/docs/sifts), and optimal growth
> temperatures from TEMPURA and from Engqvist (2018; Zenodo record 1175609). The aSAMt and AlphaFlow-MD ensembles
> were generated with the authors' released code and weights (sam2 commit ce076ed; AlphaFlow commit 0408d7c).
> All data generated in this study that are
> needed to reproduce the reported values — summary tables for every result, per-domain statistics, the
> generated aSAMt and AlphaFlow ensembles, the domain and pair lists, and the time-stamped records
> that fixed the analyses — are deposited at Zenodo (DOI: [10.5281/zenodo.XXXXXXX]). The analysis code used to produce all
> results is available at https://github.com/gaony571/framedchain-paper and archived at Zenodo (DOI:
> [10.5281/zenodo.XXXXXXX]); the method is implemented in the framedchain Python package
> (https://github.com/gaony571/framedchain, version 0.4.2, archived at Zenodo, DOI: [10.5281/zenodo.XXXXXXX]).
> Both are released under the MIT license. A file-level map from each result to its script and data file is
> given in docs/RESULTS_MAP.md of the code repository.

## Code availability sentence for the Methods (draft)

> The analysis scripts exactly as run, with the run queues, environment records and checksum records, are
> available at https://github.com/gaony571/framedchain-paper (version 1.3.1; MIT license; Zenodo DOI:
> [10.5281/zenodo.XXXXXXX]). The reusable implementation of the space, the distance and the ensemble
> statistics is the framedchain package (version 0.4.2; MIT license; Zenodo DOI: [10.5281/zenodo.XXXXXXX]).

## Software references (for the reference list)

> An H. framedchain: framed chains for protein structures and ensembles of any length [software]. Version 0.4.2.
> Zenodo; 2026. doi:[10.5281/zenodo.XXXXXXX]

> An H. framedchain-paper: analysis code for "Comparing protein conformational ensembles across sequences and lengths" [software].
> Version 1.3.1.
> Zenodo; 2026. doi:[10.5281/zenodo.XXXXXXX]

> An H. Data for "Comparing protein conformational ensembles across sequences and lengths" [dataset].
> Zenodo; 2026. doi:[10.5281/zenodo.XXXXXXX]

## Status and remaining steps

Done (2026-10-26): both repositories are on GitHub as private repositories with version tags —
`gaony571/framedchain` (tags `v0.2.0` = the version used in the runs, `v0.3.0` to `v0.4.2`) and
`gaony571/framedchain-paper` (tags `v1.0.0` = first push, `v1.1.0` = with the round44 checks of 2026-10-28, `v1.2.0` = with the round45 checks of 2026-10-28 and 10-29, `v1.3.0` = with the round46–48 checks of 2026-10-29 and 10-30, the figure and table scripts and the manuscript map). Both were renamed from framespace and framespace-paper on
2026-10-26 (GitHub redirects the old addresses). No GitHub *release* has been published yet, on purpose: Zenodo
archives a repository only for releases published after the repository is connected.

Remaining, after the manuscript is written and posted (these need the paper or the DOIs):

1. Zenodo archives a GitHub release only for a public repository that was switched on in Zenodo (log in to
   Zenodo with GitHub, then GitHub in the account settings) before the release was published. When the
   repositories are made public, switch both on and publish a GitHub release of each from its latest `main`:
   framedchain-paper `v1.3.1` (the version already in `CITATION.cff`, the README and the manuscript reference;
   `main` changed after `v1.3.0`) and framedchain `v0.4.2` (0.4.1 with a chain-break flag and documentation; no computation changed). Zenodo mints a
   DOI for each release, and `CITATION.cff` (author, affiliation, abstract, keywords, license) supplies its
   metadata. If DOIs are needed while the repositories are still private, upload `git archive` zips of the same
   tags to Zenodo by hand with reserved DOIs instead, and then do not also switch the repositories on, since each
   route makes its own record.
2. Run `tools/pack_for_zenodo.py` on the project root (Colab with Drive mounted; first with `--dry-run` to see
   the size against the Zenodo limit), check `MANIFEST_data.tsv`, and upload the data deposit. ATLAS files
   (CC BY-NC 4.0) and the structures and tables taken from them are not redistributed, nor are the MFM reference
   simulations converted to npz, which only analyses outside the paper read; the deposit lists all of these by name
   and size, with the md5 of the files that are not simulation data.
3. Fill the DOIs in the statements above (the paper title is filled in), in both `CITATION.cff` files and in the READMEs.
4. Reviewers need access at submission, and a private GitHub repository has no shareable link (only named
   collaborators can open it). Decision of 2026-10-30: both repositories are made public before submission, so
   no private link is needed and step 1 gives the two software DOIs before submission.
5. Include the time-stamped records unchanged in the data deposit, as the Methods state.
