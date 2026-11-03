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

The scripts expect the original audit layout and source artifacts; they are not
a self-contained training package. See `docs/reproducibility.md`.
