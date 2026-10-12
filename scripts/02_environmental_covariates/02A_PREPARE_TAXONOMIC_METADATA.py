#!/usr/bin/env python3
"""Select environmental metadata by the reference accession manifest."""
import argparse
from pathlib import Path
import pandas as pd
from metadata_pipeline_common import atomic_write_csv, validate_coordinates, validate_primary_key, write_manifest

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['base','out','audit-out','manifest-out']:
        p.add_argument('--'+name,required=True)
    p.add_argument('--sample-manifest',default=str(Path(__file__).resolve().parents[2]/'config/samples/taxonomic_accessions.txt'))
    p.add_argument('--id-col',default='run_accession')
    p.add_argument('--lat-col',default='lat');p.add_argument('--lon-col',default='lon')
    a=p.parse_args()
    df=pd.read_csv(a.base,low_memory=False)
    validate_primary_key(df,a.id_col,'metadata')
    validate_coordinates(df,a.lat_col,a.lon_col,'metadata')
    ids=[s.strip() for s in Path(a.sample_manifest).read_text().splitlines() if s.strip() and not s.startswith('#')]
    if not ids or len(ids)!=len(set(ids)):raise ValueError('Manifest identifiers must be nonempty and unique.')
    missing=set(ids)-set(df[a.id_col])
    if missing:raise ValueError(f'Manifest accessions absent from metadata: {sorted(missing)}')
    keep=df[a.id_col].isin(ids)
    selected=df.loc[keep].copy()
    atomic_write_csv(selected,Path(a.out))
    atomic_write_csv(pd.DataFrame({a.id_col:df[a.id_col],'included_in_reference':keep}),Path(a.audit_out))
    write_manifest(Path(a.manifest_out),Path(__file__).name,a,[Path(a.base),Path(a.sample_manifest)],Path(a.out),len(selected),len(df))
    print(f'Reference metadata: {len(selected)} profiles')
    return 0

if __name__=='__main__':raise SystemExit(main())
