#!/usr/bin/env python3
"""Apply an externally reviewed, hash-guarded metadata patch without schema changes.

Evidence JSON contains input_sha256 and changes with run_accession,
experiment_accession, sample_accession and fields: {column: {before, after}}.
Only geography/review fields listed below can change. Existing coordinates
cannot be overwritten. This does not infer or validate sampling coordinates.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

ALLOWED = {"lat", "lon", "Locality", "coordinate_review_status", "review_flags",
           "web_evidence_urls", "web_geography_check"}

def apply_patch(path, evidence):
    original = path.read_bytes()
    if hashlib.sha256(original).hexdigest() != evidence["input_sha256"]:
        raise ValueError("Input hash differs from reviewed source")
    reader = csv.DictReader(io.StringIO(original.decode("utf-8-sig")))
    fields = reader.fieldnames
    rows = list(reader)
    ids = [r["run_accession"] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate run accessions")
    lookup = {r["run_accession"]: r for r in rows}
    changes = []
    for change in evidence["changes"]:
        row = lookup[change["run_accession"]]
        for key in ("experiment_accession", "sample_accession"):
            if row[key] != change[key]:
                raise ValueError("Sample identity mismatch")
        updates = change["fields"]
        if not set(updates).issubset(ALLOWED & set(fields)):
            raise ValueError("Unapproved or new columns in patch")
        if bool("lat" in updates) != bool("lon" in updates):
            raise ValueError("Coordinate pairs must be patched together")
        if "lat" in updates:
            if row["lat"].strip() or row["lon"].strip():
                raise ValueError("Refusing to replace existing coordinates")
            lat, lon = (float(updates[k]["after"]) for k in ("lat", "lon"))
            if not (-90 <= lat <= 90 and -180 <= lon <= 180):
                raise ValueError("Invalid coordinate range")
        for col, pair in updates.items():
            if row[col] != pair["before"]:
                raise ValueError("Before-value mismatch")
            row[col] = pair["after"]
            changes.append(dict(run_accession=row["run_accession"], field=col, **pair))
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\r\n")
    writer.writeheader()
    writer.writerows(rows)
    data = output.getvalue().encode("utf-8")
    if original.startswith(b"\xef\xbb\xbf"):
        data = b"\xef\xbb\xbf" + data
    check = list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))
    before_rows = list(csv.DictReader(io.StringIO(original.decode("utf-8-sig"))))
    assert [r["run_accession"] for r in check] == ids
    allowed_cells = {(c["run_accession"], c["field"]) for c in changes}
    for a, b in zip(before_rows, check):
        assert list(a) == list(b)
        for col in fields:
            assert a[col] == b[col] or (a["run_accession"], col) in allowed_cells
    return data, dict(rows=len(rows), columns=len(fields), changed_cells=changes,
                      input_sha256=evidence["input_sha256"],
                      output_sha256=hashlib.sha256(data).hexdigest(),
                      complete_coordinate_pairs=sum(bool(r["lat"] and r["lon"]) for r in check),
                      schema_row_order_ids_and_other_cells_preserved=True)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata", type=Path, required=True)
    p.add_argument("--evidence", type=Path, required=True)
    p.add_argument("--audit", type=Path, required=True)
    a = p.parse_args()
    data, audit = apply_patch(a.metadata, json.loads(a.evidence.read_text()))
    with a.audit.open("x") as f:
        json.dump(audit, f, ensure_ascii=False, indent=2)
        f.write("\n")
    a.metadata.write_bytes(data)

if __name__ == "__main__":
    main()

