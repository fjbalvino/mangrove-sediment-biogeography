"""One-to-one sample identity contracts; no sample matching by coordinates."""
from collections import Counter
import math

PRIMARY_KEY = "run_accession"


def empty(value):
    return str(value).strip().lower() in {"", "nan", "na", "n/a", "null", "none"}


def equivalent(a, b):
    if empty(a) or empty(b):
        return empty(a) and empty(b)
    try:
        x, y = float(a), float(b)
        return math.isfinite(x) and math.isfinite(y) and math.isclose(
            x, y, rel_tol=1e-9, abs_tol=1e-8
        )
    except (TypeError, ValueError):
        return a == b


def index_by_id(rows, label, expected_rows=None):
    if expected_rows is not None and len(rows) != expected_rows:
        raise ValueError(f"{label}: expected {expected_rows} samples, found {len(rows)}")
    ids = [row.get(PRIMARY_KEY) for row in rows]
    if any(not isinstance(x, str) or empty(x) or x != x.strip() for x in ids):
        raise ValueError(f"{label}: missing, invalid or whitespace-padded {PRIMARY_KEY}")
    duplicate = sorted(x for x, n in Counter(ids).items() if n > 1)
    if duplicate:
        raise ValueError(f"{label}: duplicate {PRIMARY_KEY}: {duplicate}")
    return dict(zip(ids, rows))


def transfer_by_id(master_rows, target_rows, environmental_fields, excluded=()):
    """Left transfer of blank fields. Preserve target order, columns and row count."""
    master = index_by_id(master_rows, "master")
    target = index_by_id(target_rows, "target")
    unknown = sorted(set(target) - set(master))
    if unknown:
        raise ValueError(f"Target IDs absent from master: {unknown}")
    excluded = set(excluded)
    result = []
    for original in target_rows:
        run = original[PRIMARY_KEY]
        donor = master[run]
        row = original.copy()
        for field in ("lat", "lon"):
            if field not in original or field not in donor:
                raise ValueError(f"{run}: missing coordinate column {field}")
            if not equivalent(original[field], donor[field]):
                raise ValueError(f"{run}: coordinate conflict in {field}")
        for field in environmental_fields:
            if field not in row or field not in donor:
                raise ValueError(f"{run}: missing existing field {field}")
            if not empty(row[field]) and not empty(donor[field]):
                if not equivalent(row[field], donor[field]):
                    raise ValueError(f"{run}: conflicting existing value for {field}")
            elif run not in excluded and empty(row[field]) and not empty(donor[field]):
                row[field] = donor[field]
        result.append(row)
    if [r[PRIMARY_KEY] for r in result] != list(target):
        raise RuntimeError("ID transfer changed target membership or order")
    return result

