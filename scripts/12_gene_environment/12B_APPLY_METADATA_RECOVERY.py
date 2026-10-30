#!/usr/bin/env python3
"""Fill supported missing coordinates while preserving the original metadata schema.

Published coordinates can fill a missing pair. Figure-derived approximations
require --allow-approximate and external evidence with an uncertainty basis.
Gazetteer coordinates are contextual fields and are never promoted.
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
STATUSES = {'published_site_coordinate', 'locality_resolved_position_pending',
            'figure_georeferenced_approximate'}


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


def merge(metadata_fields, metadata, evidence_fields, evidence, allow_approximate=False):
    if not REQUIRED.issubset(evidence_fields):
        raise ValueError('Recovery table lacks required provenance fields')
    if not {'run_accession', 'experiment_accession', 'lat', 'lon'}.issubset(metadata_fields):
        raise ValueError('Metadata lacks ID or coordinate fields')
    if any(c.startswith('recovery_') for c in metadata_fields):
        raise ValueError('Use the original metadata schema, without recovery columns')
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
        if status == 'figure_georeferenced_approximate':
            if not allow_approximate:
                raise ValueError('Figure approximations require explicit --allow-approximate')
            if not row.get('recovery_uncertainty_basis') or not row.get('recovery_georeferencing_evidence'):
                raise ValueError('Figure approximations require external registration and uncertainty evidence')
            radius = numeric(row['recovery_coordinate_uncertainty_m'])
            if radius is None or radius <= 0:
                raise ValueError('Figure approximations require a positive uncertainty radius')
        if status in {'published_site_coordinate', 'figure_georeferenced_approximate'}:
            if lat is None or not (-90 <= lat <= 90 and -180 <= lon <= 180):
                raise ValueError('Recovered coordinate missing or out of range')
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
        if row['run_accession'] in changed:
            row['lat'], row['lon'] = recovered['recovery_lat'], recovered['recovery_lon']
        output.append(row)
    return metadata_fields, output, changed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--allow-approximate', action='store_true',
                        help='Allow reviewed figure-derived coordinates with external uncertainty evidence')
    args = parser.parse_args()
    evidence_fields, evidence = read_table(args.evidence)
    fields, rows, changes = merge(*read_table(args.metadata), evidence_fields, evidence,
                                 allow_approximate=args.allow_approximate)
    args.out_dir.mkdir(parents=True, exist_ok=False)
    output = args.out_dir / 'metadata_functional.csv'
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
        'original_column_order_and_schema_preserved': True,
        'provenance_kept_in_external_evidence_table': True,
        'approximate_coordinate_runs': [r['run_accession'] for r in evidence
                                       if r['recovery_status'] == 'figure_georeferenced_approximate'],
        'approximate_coordinates_explicitly_enabled': args.allow_approximate,
        'spatial_validity_uncertainty_and_raster_extraction': 'not_evaluated_by_this_script',
        'recovered_measurements_are_not_harmonized_with_spatial_proxies': True,
    }
    (args.out_dir / 'recovery_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
