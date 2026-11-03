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

## Scope

- No model was retrained.
- No historical execution log was reconstructed.
- No source PDF or third-party model file is redistributed.
- Exact sequence identity is used; no homology or subsequence substitution.
