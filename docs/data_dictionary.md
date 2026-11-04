# Data dictionary

## `data/tables/Table1_negative_evidence.csv`

One row per evidence type: missing annotation, other-function comparator,
random/generated control, measured non-response and matched food/process
control. Columns record the analytical level, what is observed, the immediate
comparison and the additional observation needed for a broader claim.

## `data/tables/Table2_reference_relationships.csv`

One row per reference-set family or unresolved row. `models` lists audit
identifiers. `sequence_verified` is the completeness status of the specified
negative version. `n_neg_available` is the number of retrieved unique
negatives; `n_neg_total` is populated only for complete versions.
`relationship` and `released_roles` are separate fields.

## `data/tables/Table3_conflict_identities.csv`

Seven sequence identities with opposite binary labels in the F665 reference.
Columns give the conflict type, H1 file positions and cross-version status.

## `data/tables/Table4_use_specific_judgments.csv`

Four intended uses. Each row states a directly interpretable comparison and
the further comparison needed to support a broader claim.

## `data/tables/Table_S1_literature_evidence.csv`

Supporting literature used for interpretation. `source_url` is a DOI link when
available. `source_file_sha256` records the examined file hash; the file itself
is not redistributed.

## `data/tables/Table_S2_audit_source_inventory.csv`

Source artifacts for the 21 study records: artifact type, source location,
version or snapshot, verification scope, byte size and SHA-256.

## `data/tables/Table_S3_audited_studies.csv`

The 21 study records, task, negative-set status, relationship, split identity
and reported negative count. `Y`, `P` and `N` are audit status codes defined
in the manuscript's Methods section.

## `data/tables/Table_S4_role_transitions.csv`

Source-to-role intersections for A2, A9 and A13, separated by label.

## `data/tables/Table_S5_conflict_versions.csv`

Version-by-sequence records for the seven conflict identities. The row count
is not a count of independent errors.

## `data/tables/Table_S5b_conflict_occurrences.csv`

Exact file and line occurrences for the seven conflict identities.

## `data/tables/Table_S6_partition_counts.csv`

Per-file record and unique-sequence counts, including known file-level
discrepancies.

## `data/tables/Table_S7_exact_overlap.csv`

Complete exact-sequence overlap matrix for the 17 retrievable negative
versions. Values are counts of identical full sequences.

## `data/checks/`

Small JSON/text outputs from the A4 loader check, the pinned Bilingua commit,
the H6 reconstruction, the software environment and the citation-key map.

## `data/evidence/Paragraph_evidence_map.csv`

Paragraph-level citation keys and evidence roles. Local source paths are
deliberately omitted.
