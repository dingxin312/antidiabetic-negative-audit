# Antidiabetic negative-evidence audit

[![DOI](https://zenodo.org/badge/1445092398.svg)](https://doi.org/10.5281/zenodo.23976752)

This repository contains the derived tables, audit code and figures supporting the manuscript:

**Negative labels in antidiabetic-peptide prediction: Benchmark reuse and implications for food-peptide screening**

The study asks how negative labels are defined, reused and evaluated in antidiabetic peptide prediction, and what a benchmark score can support in food-peptide candidate selection.

Current working version: **v0.7** (7 November 2026). The Zenodo badge and `CITATION.cff` point to the earlier v1.0.0 archive; a new archive for v0.7 is pending.

## What is included

- `data/tables/`: main and supplementary tables, including the F665 milk-precursor matches (S8), A5 candidate and tripeptide records (S9 and S9b), five candidate-selection flows (S10) and two conflict-identity assay matches (S11).
- `data/checks/`: machine-readable checks for the A4 loader, the Bilingua commit, H6 reconstruction, A5 candidates, milk-precursor counts, the software environment and citation mapping.
- `data/evidence/`: paragraph-to-citation map with citation keys and no local file paths.
- `data/reference/`: frozen UniProt FASTA and retrieval manifest for the six bovine milk proteins.
- `code/`: audit scripts retained as executable documentation. They are not a predictor package.
- `figures/`: editable PDF and PNG versions of Figures 1–3.
- `manuscript/`: v0.7 manuscript Markdown, Word file and Highlights.
- `docs/`: field definitions, reproduction notes and source policy.

## What is not included

Original journal PDFs, third-party model files and source-model repositories are not redistributed. Their identifiers, source locations and SHA-256 hashes are recorded in `data/tables/Table_S2_audit_source_inventory.csv` and in the manuscript references. Users should obtain those materials from the original sources and observe their licenses.

## Quick start

Python 3.13 or later is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Most scripts are audit snapshots and expect source artifacts at paths described in `docs/reproducibility.md`. The standalone milk-precursor check reproduces Table S8 from the four archived A3 files and the bundled frozen UniProt reference; normal reruns do not access the network. No script trains a predictor.

## Data and code availability

The current working release is this GitHub repository. Zenodo v1.0.0 (https://doi.org/10.5281/zenodo.23976753) archives the earlier v0.4 release; the v0.7 archive and version DOI are pending.

The derived tables and figures are available under the licenses in `LICENSE` and `DATA_LICENSE.md`.

## Citation

Cite the associated manuscript. The archived v1.0.0 DOI is recorded in `CITATION.cff`; update it after the v0.7 archive is created.

## Contact

Corresponding author information will be added to the manuscript before submission.
