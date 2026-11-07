# Audit code

These files are retained as executable documentation of the sequence parsing,
partition reconstruction and implementation checks described in the
manuscript.

- `parse_negative_sets.py`: parses released sequence files, normalizes records
  and writes long-form sequence tables.
- `reconstruct_and_inspect.py`: reconstructs H6 80/20 splits for seeds 1–10 and
  executes isolated A4 loader functions.
- `independent_recheck_claims.py`: recomputes selected membership and conflict
  claims from parsed records.
- `validate_materials.py`: checks counts, partitions, conflict identities,
  source hashes and table assembly.
- `milk_precursor_audit.py`: checks the four archived A3 F665 input hashes and
  the frozen UniProt reference hash, then reproduces Table S8 and its summary.
- `a5_candidate_audit.py`: extracts the A5 pentapeptide and tripeptide tables from the archived ACS supporting PDF and checks their F665 membership. It is an audit snapshot, not a standalone predictor package.

The scripts expect the original audit layout and source artifacts; they are not
a self-contained training package. See `docs/reproducibility.md`.
The milk-precursor script is standalone once the four A3 files are obtained.
