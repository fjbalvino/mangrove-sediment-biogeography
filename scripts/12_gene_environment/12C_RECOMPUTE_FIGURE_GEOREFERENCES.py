#!/usr/bin/env python3
"""Recompute map-marker coordinates from external, visually reviewed image matches.

Requires NumPy. Does not infer field GPS positions, locate mangrove pixels,
extract environmental rasters or modify metadata. Input contains panels with
source_points/reference_points and the actual ArcGIS EPSG:3857 export extent,
plus sites with panel/source_center_px/source_marker_radius_px. All image pixels
are zero-based; 0.5 offsets map pixel centers to the reference extent.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np


def recompute(evidence):
    diagnostics, sites = {}, []
    for key, panel in evidence['panels'].items():
        p = np.asarray(panel['source_points'], dtype=float)
        q = np.asarray(panel['reference_points'], dtype=float)
        if p.shape != q.shape or p.ndim != 2 or p.shape[1] != 2 or not np.isfinite([p,q]).all():
            raise ValueError('Invalid paired image matches')
        _, indices = np.unique(np.round(p, 2), axis=0, return_index=True)
        p, q = p[indices], q[indices]
        A = np.c_[p, np.ones(len(p))]
        if len(p) < 6 or np.linalg.matrix_rank(A) != 3:
            raise ValueError('At least six non-collinear unique image matches required')
        M = np.linalg.lstsq(A, q, rcond=None)[0].T
        residual = np.linalg.norm(A @ M.T - q, axis=1)
        loo = []
        for i in range(len(p)):
            keep = np.arange(len(p)) != i
            if np.linalg.matrix_rank(A[keep]) != 3:
                raise ValueError('Degenerate leave-one-out control set')
            fit = np.linalg.lstsq(A[keep], q[keep], rcond=None)[0]
            loo.append(float(np.linalg.norm(A[i] @ fit - q[i])))
        export = panel['reference_export']
        e = export['extent']
        sr = e['spatialReference']
        if sr.get('latestWkid', sr.get('wkid')) not in {3857, 102100}:
            raise ValueError('Reference must be Web Mercator EPSG:3857')
        dx = (e['xmax']-e['xmin']) / export['width']
        dy = (e['ymax']-e['ymin']) / export['height']
        if dx <= 0 or dy <= 0:
            raise ValueError('Invalid reference extent or dimensions')
        diagnostics[key] = dict(matrix=M.tolist(), unique_matched_features=len(p),
                                rmse_pixels=float(np.sqrt(np.mean(residual**2))),
                                loo_max_pixels=max(loo),
                                validation='Internal registration diagnostics, not field-position accuracy')
        for site in (s for s in evidence['sites'] if s['panel'] == key):
            px, py = M @ np.r_[site['source_center_px'], 1]
            X = e['xmin'] + (px+.5)*dx
            Y = e['ymax'] - (py+.5)*dy
            lon = math.degrees(X/6378137)
            lat = math.degrees(2*math.atan(math.exp(Y/6378137))-math.pi/2)
            ground = math.cos(math.radians(lat))
            radius = float(site['source_marker_radius_px'] *
                           np.linalg.svd(np.diag([dx,dy]) @ M[:,:2], compute_uv=False).max() * ground)
            sites.append(dict(label=site['label'], panel=key, lat=round(lat,3), lon=round(lon,3),
                              lat_unrounded=lat, lon_unrounded=lon,
                              marker_max_radius_ground_m=radius,
                              status='figure_georeferenced_approximate'))
    return dict(panels=diagnostics, sites=sites,
                warning='Original author marker placement and composite-sample extent remain unknown.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = recompute(json.loads(args.evidence.read_text()))
    result['input_sha256'] = hashlib.sha256(args.evidence.read_bytes()).hexdigest()
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2, ensure_ascii=False)
        handle.write('\n')


if __name__ == '__main__':
    main()
