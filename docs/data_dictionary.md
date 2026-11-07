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
Columns give the conflict type, H1 FASTA source-line positions, cross-version
status and one-based positions in the P02666 precursor. Source-line positions
are not record ordinals; a P02666 substring is not evidence of actual release.

## `data/tables/Table4_use_specific_judgments.csv`

Four intended uses. Columns separate the audited sequence fact, the direct
model or food observation, and the next comparison for a broader claim.

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

## `data/tables/Table_S8_milk_precursor_matches.csv`

One row per sequence identity within each F665 label class (665 positive and
665 negative rows). The seven identities in both classes appear twice, once
per class. Columns record split, length, exact matches to six selected bovine
milk precursors, one-based P02666 start positions and opposite-label status.
No column asserts actual milk origin, digestion release, activity or model use.

`data/reference/Reference_milk_proteins_20261105.fasta` freezes six reviewed
UniProt sequences; its JSON retrieval manifest lists accession, sequence
version, source URL, time and SHA-256. `data/checks/milk_precursor_summary.json`
contains length-stratified counts, source hashes, seven conflict positions and
the >=5-residue main and >=6-residue sensitivity calculations.

## `data/checks/`

Small JSON/text outputs from the A4 loader check, the pinned Bilingua commit,
the H6 reconstruction, milk precursor matches, the software environment and
the citation-key map.

## `data/evidence/Paragraph_evidence_map.csv`

Paragraph-level citation keys and evidence roles. Local source paths are
deliberately omitted.


## `data/tables/Table_S9_A5_candidate_results.csv`

A5 five-residue design: sequence, predicted class, reported IC50 and source page. Eight designs had a negative prediction and a measured IC50; IPIPI is the lowest measured IC50 in the set.

## `data/tables/Table_S9b_A5_tripeptide_membership.csv`

A5 tripeptide designs with predicted class and exact membership in the F665 positive or negative files. The 80 tripeptides were all predicted positive; 72 had a measured IC50.

## `data/tables/Table_S10_candidate_selection_stages.csv`

Five purposefully selected food-protein candidate flows. Columns separate the enumerated library, the filtered set and the assayed peptides. Rows are not pooled into a recovery or false-negative rate.

## `data/tables/Table_S11_conflict_assay_matches.csv`

Two exact F665 conflict identities with independently published in vitro DPP-IV assay values. The assay result does not establish the origin of the deposited negative label.
