"""Rebuild A5 candidate counts from the archived original ACS supporting PDF."""

import csv
import hashlib
import json
import re
import subprocess
from collections import Counter
from pathlib import Path


SUPPLEMENT = Path(__file__).resolve().parents[1]
AUDIT_ROOT = SUPPLEMENT.parents[1]
SOURCE = AUDIT_ROOT / "A5_BERT-DPPIV/source/ao3c05571_si_001.pdf"
EXPECTED_SHA256 = "99407648164cf92e2828a381ceff62e7723b60dec8e09f615ab8c57b8bc762f4"


def page(number):
    return subprocess.check_output(
        ["pdftotext", "-f", str(number), "-l", str(number), "-layout", str(SOURCE), "-"],
        text=True,
    )


assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == EXPECTED_SHA256
tripeptides = page(8)
assert "Supplementary Table. 3" in tripeptides
tri_rows = re.findall(
    r"^\s*([A-Z]{3})\s+([01])\s+([A-Z]{3})\s+([01])"
    r"\s+([A-Z]{3})\s+([01])\s+([A-Z]{3})\s+([01])\s*$",
    tripeptides,
    re.M,
)
tri = {name: int(prediction) for row in tri_rows for name, prediction in zip(row[::2], row[1::2])}
assert len(tri_rows) == 20 and len(tri) == 80 and set(tri.values()) == {1}

prediction_text = page(9)
assert "Supplementary Table. 4" in prediction_text
prediction_rows = re.findall(
    r"^\s*([A-Z]{5})\s+([01])\s+([A-Z]{5})\s+([01])\s*$",
    prediction_text,
    re.M,
)
predictions = {
    name: int(prediction)
    for row in prediction_rows
    for name, prediction in zip(row[::2], row[1::2])
}
assert len(prediction_rows) == 11 and len(predictions) == 22

with (SUPPLEMENT / "Table_S8_milk_precursor_matches.csv").open(encoding="utf-8", newline="") as handle:
    f665_rows = list(csv.DictReader(handle))
assert len(f665_rows) == 1330
f665_roles = {}
for row in f665_rows:
    f665_roles.setdefault(row["sequence"], []).append(f"{row['label']}:{row['split']}")
tri_overlap = Counter(role for sequence in tri for role in f665_roles.get(sequence, ()))
assert tri_overlap == {"positive:train": 12, "positive:test": 6}
assert all(sequence not in f665_roles for sequence in predictions)

measurement_text = page(10)
assert "Supplementary Table. 5" in measurement_text
measurement_rows = re.findall(
    r"^\s*([A-Z]{5})\s+(\d+(?:\.\d+)?)\s+-\d+(?:\.\d+)?"
    r"\s+([A-Z]{5})\s+(\d+(?:\.\d+)?)\s+-\d+(?:\.\d+)?\s*$",
    measurement_text,
    re.M,
)
ic50 = {name: value for row in measurement_rows for name, value in zip(row[::2], row[1::2])}
assert len(measurement_rows) == 11 and predictions.keys() == ic50.keys()
assert sum(prediction == 0 for prediction in predictions.values()) == 8
assert predictions["IPIPI"] == 0 and float(ic50["IPIPI"]) == 47.47
assert min(ic50, key=lambda peptide: float(ic50[peptide])) == "IPIPI"

output = SUPPLEMENT / "Table_S9_A5_candidate_results.csv"
with output.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["sequence", "SI_Table_S4_prediction", "SI_Table_S5_IC50_uM", "prediction_page", "measurement_page"])
    for peptide in sorted(predictions):
        writer.writerow([peptide, predictions[peptide], ic50[peptide], 9, 10])

with (SUPPLEMENT / "Table_S9b_A5_tripeptide_membership.csv").open(
    "w", encoding="utf-8", newline=""
) as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["sequence", "SI_Table_S3_prediction", "F665_audited_version_label_and_split", "prediction_page"])
    for peptide in sorted(tri):
        writer.writerow([peptide, tri[peptide], ";".join(f665_roles.get(peptide, ())), 8])

summary = {
    "source_relative_to_audit_root": "A5_BERT-DPPIV/source/ao3c05571_si_001.pdf",
    "source_sha256": EXPECTED_SHA256,
    "f665_reference_sha256": hashlib.sha256(
        (SUPPLEMENT / "Table_S8_milk_precursor_matches.csv").read_bytes()
    ).hexdigest(),
    "tripeptides_predicted_positive_in_S3": len(tri),
    "tripeptides_exact_in_F665": {"total": sum(tri_overlap.values()), **dict(sorted(tri_overlap.items()))},
    "pentapeptides_in_S4_and_S5": len(predictions),
    "pentapeptides_exact_in_F665": 0,
    "pentapeptides_predicted_negative_with_reported_IC50": sum(value == 0 for value in predictions.values()),
    "lowest_reported_IC50_uM": {"sequence": "IPIPI", "value": 47.47, "predicted": 0},
    "limits": "Exact membership uses the audited A3 F665 positive and negative version (S8), with A5 sequence-file identity checked separately in the source audit. S3-S5 record these designed libraries; 72/80 measurable tripeptide responses appear in article Table 1. Sequence overlap does not resolve assay conditions, and zero pentapeptide overlap is not a homology test. No F665 error attribution or external-library error rate is estimated.",
}
(SUPPLEMENT / "a5_candidate_summary.json").write_text(
    json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
)
print(f"Rebuilt {len(tri)} tripeptide predictions (18 F665 matches: 12 train, 6 test) and {len(predictions)} pentapeptide rows (0 F665 matches); 8 predicted negatives; IPIPI 47.47 uM.")
