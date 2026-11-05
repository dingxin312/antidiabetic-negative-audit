"""Reproduce the bounded F665 length and bovine-milk precursor annotation.

Run with --fetch-reference only to prepare the dated UniProt FASTA snapshot.
Ordinary reruns use that frozen FASTA and never contact the network.
"""

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from urllib.request import urlopen


ACCESSIONS = ("P02662", "P02663", "P02666", "P02668", "P02754", "P00711")
SOURCE_SHA256 = {
    "TR_pos_dppiv.txt": "669d3ba252705af989b44ef14c1ea3429870e21e00b3f096fb73fd779aee558e",
    "TS_pos_dppiv.txt": "e96696b631df341976fad2b688c681d66976cb8a7229ceb5d9dc418a20ffd065",
    "TR_neg_dppiv.txt": "5ab5a88dca0bdde26d835240d6923eac4a7fe239e2ed4c35c5904cd1a649d113",
    "TS_neg_dppiv.txt": "7094b9386402a5e566fe39fb3f73120d32753ccf76cfc91a9f12918e2eb97444",
}
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
LENGTH_BINS = ((2, 4), (5, 5), (6, 7), (8, 9), (10, 14), (15, 19), (20, None))


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def parse_fasta(text):
    header, parts = None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(parts).upper()
            header, parts = line, []
        elif header is None:
            raise ValueError("Sequence before FASTA header")
        else:
            parts.append(line)
    if header is not None:
        yield header, "".join(parts).upper()


def checked_sequence(sequence, source):
    if len(sequence) < 2 or set(sequence) - STANDARD_AA:
        raise ValueError(f"Non-standard peptide/protein sequence in {source}: {sequence!r}")
    return sequence


def fetch_reference(path):
    if path.exists():
        raise FileExistsError(f"Refusing to replace frozen reference: {path}")
    blocks, records = [], []
    for accession in ACCESSIONS:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
        with urlopen(url, timeout=30) as response:
            block = response.read()
        entries = list(parse_fasta(block.decode("utf-8")))
        if len(entries) != 1 or not entries[0][0].startswith(f">sp|{accession}|"):
            raise ValueError(f"Unexpected UniProt record for {accession}")
        checked_sequence(entries[0][1], accession)
        blocks.append(block.rstrip(b"\n") + b"\n")
        records.append({"accession": accession, "url": url, "header": entries[0][0],
                        "fasta_sha256": sha256(block), "sequence_length": len(entries[0][1])})
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"".join(blocks))
    metadata = {"retrieved_utc": datetime.now(timezone.utc).isoformat(),
                "reference_file": path.name, "reference_sha256": sha256(path.read_bytes()),
                "records": records}
    path.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def read_sources(source_dir):
    labels = defaultdict(dict)
    for name, expected in SOURCE_SHA256.items():
        path = source_dir / name
        raw = path.read_bytes()
        if sha256(raw) != expected:
            raise ValueError(f"A3 source hash differs from audited snapshot: {path}")
        label = "positive" if "_pos_" in name else "negative"
        split = "train" if name.startswith("TR_") else "test"
        for header, sequence in parse_fasta(raw.decode("utf-8-sig")):
            checked_sequence(sequence, name)
            if not re.fullmatch(fr">{label.title()} \d+", header):
                raise ValueError(f"Unexpected source FASTA header in {name}: {header}")
            if sequence in labels[label]:
                raise ValueError(f"Repeated identity within {label}: {sequence}")
            labels[label][sequence] = split
    for label in ("positive", "negative"):
        expected = {"train": 532, "test": 133}
        observed = {split: sum(role == split for role in labels[label].values()) for split in expected}
        if observed != expected:
            raise ValueError(f"Expected F665 {label} train/test counts {expected}, got {observed}")
    return labels


def positions(sequence, protein):
    found, start = [], 0
    while (start := protein.find(sequence, start)) >= 0:
        found.append(start + 1)
        start += 1
    return found


def run(source_dir, reference_path, output_dir):
    reference_bytes = reference_path.read_bytes()
    reference = {}
    for header, sequence in parse_fasta(reference_bytes.decode("utf-8")):
        match = re.match(r">sp\|([^|]+)\|", header)
        if not match or match.group(1) in reference:
            raise ValueError(f"Unexpected reference FASTA header: {header}")
        reference[match.group(1)] = checked_sequence(sequence, header)
    if tuple(reference) != ACCESSIONS:
        raise ValueError(f"Expected six ordered, reviewed UniProt entries: {ACCESSIONS}")
    metadata_path = reference_path.with_suffix(".json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata["reference_sha256"] != sha256(reference_bytes):
        raise ValueError("Frozen UniProt reference hash differs from its retrieval manifest")
    for entry in metadata["records"]:
        accession = entry["accession"]
        if (accession not in reference or len(reference[accession]) != entry["sequence_length"]
                or not entry["header"].startswith(f">sp|{accession}|")):
            raise ValueError(f"Frozen UniProt entry differs from its retrieval manifest: {accession}")

    labels = read_sources(source_dir)
    conflicts = labels["positive"].keys() & labels["negative"].keys()
    expected_conflicts = {"FAQTQS", "GPFPIIV", "LHLPLP", "LPVPQ", "LTDVEN", "LTQTPVVVPPF", "VRGPFP"}
    if conflicts != expected_conflicts:
        raise ValueError(f"Unexpected F665 opposite-label identities: {sorted(conflicts)}")

    rows, summary = [], {"classes": {}, "length_strata": {}, "conflicts": []}
    for label in ("positive", "negative"):
        for sequence, split in sorted(labels[label].items()):
            match_positions = {acc: positions(sequence, protein) for acc, protein in reference.items()}
            matches = {acc: loc for acc, loc in match_positions.items() if loc}
            rows.append({"label": label, "split": split, "sequence": sequence,
                         "length": len(sequence), "matched_accessions": ";".join(matches),
                         "beta_casein_start_1based": ";".join(map(str, matches.get("P02666", []))),
                         "opposite_label_identity": int(sequence in conflicts)})
        subset = [r for r in rows if r["label"] == label]
        ge5 = [r for r in subset if r["length"] >= 5]
        ge6 = [r for r in subset if r["length"] >= 6]
        summary["classes"][label] = {
            "identities": len(subset), "median_length": median(r["length"] for r in subset),
            "length_at_most_4": sum(r["length"] <= 4 for r in subset),
            "length_at_least_5": len(ge5),
            "six_protein_matches_length_at_least_5": sum(bool(r["matched_accessions"]) for r in ge5),
            "beta_casein_matches_length_at_least_5": sum(bool(r["beta_casein_start_1based"]) for r in ge5),
            "length_at_least_6": len(ge6),
            "six_protein_matches_length_at_least_6": sum(bool(r["matched_accessions"]) for r in ge6),
        }
        summary["length_strata"][label] = [
            {"length": f"{low}+" if high is None else f"{low}-{high}",
             "identities": len(group := [r for r in subset if r["length"] >= low and (high is None or r["length"] <= high)]),
             "six_protein_matches": sum(bool(r["matched_accessions"]) for r in group),
             "beta_casein_matches": sum(bool(r["beta_casein_start_1based"]) for r in group)}
            for low, high in LENGTH_BINS
        ]
    for sequence in sorted(conflicts):
        matches = positions(sequence, reference["P02666"])
        if not matches:
            raise ValueError(f"Conflict identity not in P02666: {sequence}")
        summary["conflicts"].append({"sequence": sequence,
                                     "beta_casein_precursor_start_1based": matches,
                                     "positive_split": labels["positive"][sequence],
                                     "negative_split": labels["negative"][sequence]})
    summary.update({"scope": "A3 archived copy of F665; exact full-sequence matches to six selected Bos taurus reference precursors",
                    "source_files_sha256": SOURCE_SHA256,
                    "reference_file": reference_path.name,
                    "reference_sha256": sha256(reference_bytes),
                    "threshold": "Matches presented for length >=5; >=6 sensitivity and all-length bins included",
                    "cross_label_overlap": len(conflicts),
                    "interpretation": "Sequence containment is compatible with a precursor, not evidence of actual peptide origin, release, activity or model reliance on precursor identity."})
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "Table_S8_milk_precursor_matches.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (output_dir / "milk_precursor_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"classes": summary["classes"], "conflicts": summary["conflicts"]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", type=Path, help="Path to A3_StackDPPIV/source containing the four TR/TS FASTA files")
    parser.add_argument("--reference", type=Path, default=Path(__file__).resolve().parents[1] / "data" / "reference" / "Reference_milk_proteins_20261105.fasta")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--fetch-reference", action="store_true", help="Fetch and freeze six reviewed UniProt FASTA entries once")
    args = parser.parse_args()
    if args.fetch_reference:
        fetch_reference(args.reference)
    if not args.fetch_reference:
        run(args.source_dir, args.reference, args.output)
