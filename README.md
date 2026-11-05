# Antidiabetic negative-evidence audit

[![DOI](https://zenodo.org/badge/1445092398.svg)](https://doi.org/10.5281/zenodo.23976752)

This repository contains the derived tables, audit code and figures supporting the manuscript:

**Negative evidence in AI-assisted discovery of food-derived antidiabetic peptides: construction, benchmark reuse and validation boundaries**

The study asks how negative labels are defined, inherited and evaluated in antidiabetic peptide prediction, and what a benchmark score can support in food-peptide screening.

## What is included

- `data/tables/`: derived tables used in the manuscript and supplementary material.
- `data/checks/`: small machine-readable checks for the A4 loader, the Bilingua commit, H6 reconstructions, the software environment and citation mapping.
- `data/evidence/`: paragraph-to-citation map with source keys, without local file paths.
- `code/`: audit scripts retained as executable documentation. They are not a predictor package.
- `figures/`: editable PDF and PNG versions of Figures 1–3.
- `manuscript/`: manuscript Markdown, Word file and Highlights.
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

The audit code and derived tables are available at https://doi.org/10.5281/zenodo.23976753.

The derived tables and figures are available under the licenses in `LICENSE` and `DATA_LICENSE.md`.

## Citation

Cite the associated manuscript. The archived DOI is recorded in `CITATION.cff`.

## Contact

Corresponding author information will be added to the manuscript before submission.
