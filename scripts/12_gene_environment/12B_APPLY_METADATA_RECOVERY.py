#!/usr/bin/env python3
"""Join externally curated literature evidence without replacing measured/proxy fields.

Only `recovery_status=published_site_coordinate` can fill a missing coordinate
pair. Gazetteer coordinates are contextual fields and are never promoted.
No extraction, geocoding, spatial validation or unit harmonization is performed.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

REQUIRED = {
    'run_accession', 'experiment_accession', 'recovery_status',
    'recovery_source_url', 'recovery_identity_basis',
    'recovery_lat', 'recovery_lon', 'recovery_locality',
    'recovery_coordinate_uncertainty_m', 'recovery_datum',
}
STATUSES = {'published_site_coordinate', 'locality_resolved_position_pending'}


def read_table(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        if not fields or len(fields) != len(set(fields)):
            raise ValueError('Missing or duplicate column names')
        rows = list(reader)
    if any(None in row or None in row.values() for row in rows):
        raise ValueError('Malformed CSV row')
    ids = [row['run_accession'] for row in rows]
    if not all(ids) or len(ids) != len(set(ids)):
        raise ValueError('Missing or duplicate run_accession')
    return fields, rows


def numeric(value):
    if value.strip().lower() in {'', 'na', 'nan', 'none', 'null'}:
        return None
    number = float(value)
    if not math.isfinite(number):
        raise ValueError('Non-finite coordinate or uncertainty')
    return number


def merge(metadata_fields, metadata, evidence_fields, evidence):
    if not REQUIRED.issubset(evidence_fields):
        raise ValueError('Recovery table lacks required provenance fields')
    if not {'run_accession', 'experiment_accession', 'lat', 'lon'}.issubset(metadata_fields):
        raise ValueError('Metadata lacks ID or coordinate fields')
    added = [c for c in evidence_fields if c.startswith('recovery_')]
    if set(added) & set(metadata_fields):
        raise ValueError('Recovery fields already exist; refusing implicit replacement')
    lookup = {r['run_accession']: r for r in metadata}
    changed = []
    for row in evidence:
        key = row['run_accession']
        if key not in lookup or row['experiment_accession'] != lookup[key]['experiment_accession']:
            raise ValueError('Recovery identity mismatch: ' + key)
        status = row['recovery_status']
        if status not in STATUSES or not row['recovery_source_url'] or not row['recovery_identity_basis']:
            raise ValueError('Missing source/identity basis or unsupported recovery status')
        lat, lon = numeric(row['recovery_lat']), numeric(row['recovery_lon'])
        if (lat is None) != (lon is None):
            raise ValueError('Incomplete recovered coordinate pair')
        if status == 'published_site_coordinate':
            if lat is None or not (-90 <= lat <= 90 and -180 <= lon <= 180):
                raise ValueError('Published coordinate missing or out of range')
            if not row['recovery_datum'] or not row['recovery_locality']:
                raise ValueError('Datum status and locality are required')
            old = lookup[key]
            if numeric(old['lat']) is not None or numeric(old['lon']) is not None:
                raise ValueError('Existing coordinates cannot be overwritten: ' + key)
            uncertainty = numeric(row['recovery_coordinate_uncertainty_m'])
            if uncertainty is not None and uncertainty <= 0:
                raise ValueError('Unknown uncertainty must be blank, not zero')
            changed.append(key)
        elif lat is not None:
            raise ValueError('Unresolved position cannot supply sampling coordinates')
    by_id = {r['run_accession']: r for r in evidence}
    output = []
    for original in metadata:
        row = dict(original)
        recovered = by_id.get(row['run_accession'])
        row.update({c: recovered[c] if recovered else '' for c in added})
        if row['run_accession'] in changed:
            row['lat'], row['lon'] = recovered['recovery_lat'], recovered['recovery_lon']
        output.append(row)
    return metadata_fields + added, output, changed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    fields, rows, changes = merge(*read_table(args.metadata), *read_table(args.evidence))
    args.out_dir.mkdir(parents=True, exist_ok=False)
    output = args.out_dir / 'metadata_191_with_literature_recovery.csv'
    with output.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    audit = {
        'metadata_rows': len(rows), 'recovered_coordinate_runs': changes,
        'recovered_coordinate_count': len(changes),
        'source_metadata_sha256': hashlib.sha256(args.metadata.read_bytes()).hexdigest(),
        'source_evidence_sha256': hashlib.sha256(args.evidence.read_bytes()).hexdigest(),
        'output_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
        'only_original_columns_changed': ['lat', 'lon'] if changes else [],
        'original_site_ids_localities_regions_and_basis_preserved': True,
        'use_recovery_fields_for_new_provenance': True,
        'spatial_validity_uncertainty_and_raster_extraction': 'not_evaluated_by_this_script',
        'recovered_measurements_are_not_harmonized_with_spatial_proxies': True,
    }
    (args.out_dir / 'recovery_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
