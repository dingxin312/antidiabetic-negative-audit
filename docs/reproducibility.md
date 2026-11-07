# Reproduction notes

These scripts are audit snapshots, not a self-contained predictor package.
They document how sequence files were parsed, how partitions were compared and
how the implementation checks were performed.

## Environment

- Python 3.13.9
- numpy 2.3.5
- scikit-learn 1.7.2
- openpyxl 3.1.5
- python-docx 1.1 or later

The same information is in `data/checks/environment.json`.

## Expected source layout

The scripts were written against the original audit layout. Source artifacts
are not bundled because of size and third-party rights. The inventory in
`data/tables/Table_S2_audit_source_inventory.csv` gives each artifact's source
location and SHA-256.

`reconstruct_and_inspect.py` accepts a project root:

```bash
python code/reconstruct_and_inspect.py /path/to/audit_project
```

The other scripts expect the project's `02_解析/`, `03_Table2/`, `06_独立复核/`
and `07_封版/` directories. They are included as executable documentation of
the checks reported in the manuscript, not as a turnkey pipeline.

## A5 candidate audit

`code/a5_candidate_audit.py` is retained as an executable snapshot of the A5 candidate and membership extraction. It requires the original `A5_BERT-DPPIV/source/ao3c05571_si_001.pdf` and the project layout used for the audit; the repository does not redistribute that PDF. The derived tables are `data/tables/Table_S9_A5_candidate_results.csv`, `Table_S9b_A5_tripeptide_membership.csv` and `data/checks/a5_candidate_summary.json`.

## Scope

- No model was retrained.
- No historical execution log was reconstructed.
- No source PDF or third-party model file is redistributed.
- Exact sequence identity is used; no homology or subsequence substitution.

## Frozen milk-precursor annotation

`data/reference/Reference_milk_proteins_20261105.fasta` contains the exact
UniProt snapshot used for the six bovine milk proteins. They were selected to cover the major bovine casein and whey fractions: alpha-s1-casein, alpha-s2-casein, beta-casein, kappa-casein, beta-lactoglobulin and alpha-lactalbumin (P02662,
P02663, P02666, P02668, P02754 and P00711). Its adjacent JSON file records
the retrieval time, individual reviewed FASTA headers and SHA-256 hashes.
The input is the four A3 TR/TS files listed in Table S2. Their expected hashes
are in the script and are checked before parsing. With these third-party
inputs obtained in the `A3_StackDPPIV/source` layout, reproduce offline:

```bash
python code/milk_precursor_audit.py /path/to/A3_StackDPPIV/source \
  --reference data/reference/Reference_milk_proteins_20261105.fasta \
  --output /tmp/milk_precursor_recheck
cmp /tmp/milk_precursor_recheck/Table_S8_milk_precursor_matches.csv \
  data/tables/Table_S8_milk_precursor_matches.csv
cmp /tmp/milk_precursor_recheck/milk_precursor_summary.json \
  data/checks/milk_precursor_summary.json
```

The analysis compares distinct exact sequence identities within each label
class to six selected reference precursors, with >=5 residues for the main
comparison, >=6 as a sensitivity threshold and seven length strata in the
summary. The positive and negative sets each contain 665 identities and
overlap by seven, so class denominators must not be added as if disjoint.
Substring matches do not establish experimental origin, peptide release or
predictor behavior. The one-time `--fetch-reference` option creates a fresh
FASTA and retrieval manifest only at a new path; use the included frozen file
to reproduce the reported results.
