# Availability statements and software references (PLOS Computational Biology)

## What the journal requires (code availability policy, checked 2026-10-26)

Source: https://journals.plos.org/ploscompbiol/s/code-availability

| Requirement | How it is met |
|---|---|
| Code that directly relates to the results is shared upon publication; editors and reviewers can access it at submission (a private link is acceptable) | This repository (all analysis scripts as run) and the framedchain package; both public before submission (decision of 2026-10-30), so editors and reviewers can open them at submission |
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
> that fixed the analyses — are deposited at Zenodo (DOI: 10.5281/zenodo.23858418). The analysis code used to produce all
> results is available at https://github.com/gaony571/framedchain-paper and archived at Zenodo (DOI:
> [10.5281/zenodo.XXXXXXX]); the method is implemented in the framedchain Python package
> (https://github.com/gaony571/framedchain, version 0.4.3, archived at Zenodo, DOI: [10.5281/zenodo.XXXXXXX]).
> Both are released under the MIT license. A file-level map from each result to its script and data file is
> given in docs/RESULTS_MAP.md of the code repository.

## Code availability sentence for the Methods (draft)

> The analysis scripts exactly as run, with the run queues, environment records and checksum records, are
> available at https://github.com/gaony571/framedchain-paper (version 1.3.1; MIT license; Zenodo DOI:
> [10.5281/zenodo.XXXXXXX]). The reusable implementation of the space, the distance and the ensemble
> statistics is the framedchain package (version 0.4.3; MIT license; Zenodo DOI: [10.5281/zenodo.XXXXXXX]).

## Software references (for the reference list)

> An H. framedchain: framed chains for protein structures and ensembles of any length [software]. Version 0.4.3.
> Zenodo; 2026. doi:[10.5281/zenodo.XXXXXXX]

> An H. framedchain-paper: analysis code for "Comparing protein conformational ensembles across sequences and lengths" [software].
> Version 1.3.1.
> Zenodo; 2026. doi:[10.5281/zenodo.XXXXXXX]

> An H. Data for "Comparing protein conformational ensembles across sequences and lengths" [dataset].
> Zenodo; 2026. doi:10.5281/zenodo.23858418

## Status (as of 2026-11-02, before the releases)

Both repositories are on GitHub with version tags (from 2026-10-26) —
`gaony571/framedchain` (tags `v0.2.0` = the version used in the runs, `v0.3.0` to `v0.4.2`; `v0.4.3` is created with
its release) and
`gaony571/framedchain-paper` (tags `v1.0.0` = first push, `v1.1.0` = with the round44 checks of 2026-10-28, `v1.2.0` = with the round45 checks of 2026-10-28 and 10-29, `v1.3.0` = with the round46–48 checks of 2026-10-29 and 10-30, the figure and table scripts and the manuscript map; `v1.3.1` is created with its release). Both were renamed from framespace and framespace-paper on
2026-10-26 (GitHub redirects the old addresses). They were kept private until the releases below. No GitHub *release*
had been published before them, on purpose: Zenodo archives a repository only for releases published after the
repository is connected.

Steps before submission:

1. Zenodo archives a GitHub release only for a public repository that was switched on in Zenodo (log in to
   Zenodo with GitHub, then GitHub in the account settings) before the release was published. When the
   repositories are made public, switch both on and publish a GitHub release of each:
   framedchain-paper `v1.3.1` from its latest `main` (the version already in `CITATION.cff`, the README and the
   manuscript reference; `main` changed after `v1.3.0`) and framedchain `v0.4.3` (0.4.1 with the
   chain-break flag of 0.4.2 and documentation; no computation changed). Zenodo mints a DOI for each release, and `CITATION.cff` (author, affiliation, abstract, keywords, license) supplies its
   metadata.
2. The data deposit was packed with `tools/pack_for_zenodo.py` on the project root (Colab with Drive mounted), with
   the files listed in its `MANIFEST_data.tsv` (built with `--ensembles`, so the generated aSAMt and AlphaFlow
   ensembles are included), and uploaded on 2026-11-01 (18:36 UTC; Zenodo's md5 matched for all three files). Before
   it is published, `EXCLUDED_public.tsv`, rebuilt with `--excluded-only` to add inputs the first list missed
   (`docs/DATA.md`), replaces the uploaded one, and `prespecification_document_versions.zip` (the three surviving
   versions of the pre-specification document) is added. ATLAS files
   (CC BY-NC 4.0) and the structures and tables taken from them are not redistributed, nor are the MFM reference
   simulations converted to npz, which only analyses outside the paper read; the deposit lists all of these by name
   and size, with the md5 of the files of up to 64 MiB that are not simulation data. It holds the time-stamped
   records unchanged, as the Methods state. The deposit is the Zenodo draft 23858418, with the reserved DOI
   10.5281/zenodo.23858418 (filled in above); it is published before submission.
3. Fill the DOIs in the statements above (the paper title is filled in), in both `CITATION.cff` files and in the
   READMEs. The data DOI is filled in; the two software DOIs follow the releases of step 1.
4. Reviewers need access at submission, and a private GitHub repository has no shareable link (only named
   collaborators can open it). Decision of 2026-10-30: both repositories are made public before submission, so
   no private link is needed and step 1 gives the two software DOIs before submission.
