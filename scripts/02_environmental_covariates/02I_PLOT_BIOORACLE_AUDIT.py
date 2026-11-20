#!/usr/bin/env python3
"""Render accession/coordinate coverage from a 02H or 02J salinity audit."""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch
import numpy as np

COLORS = {"native": "#187b83", "coastal": "#73bdad", "cardinal": "#496aa1", "missing": "#c77847"}


def group(status):
    if status in {"assigned_native", "assigned_reviewed_native", "assigned_recovered_native"}:
        return "native"
    if status == "assigned_reviewed_coastal":
        return "coastal"
    if status == "assigned_cardinal":
        return "cardinal"
    return "missing"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True)
    parser.add_argument("--land-geojson", required=True, help="External Natural Earth land polygons, WGS84")
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()
    with Path(args.audit).open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    counts = collections.Counter(group(r["status"]) for r in rows)
    countries = collections.defaultdict(collections.Counter)
    sites = {}
    for r in rows:
        countries[r["Country"]][group(r["status"])] += 1
        try:
            point = float(r["lat"]), float(r["lon"])
            if all(math.isfinite(v) for v in point):
                sites.setdefault(point, set()).add(group(r["status"]))
        except ValueError:
            pass
    assigned = counts["native"] + counts["coastal"] + counts["cardinal"]
    selected_sites = {key for key, states in sites.items() if states & {"native", "coastal", "cardinal"}}
    cells = {(r["cell_lat"], r["cell_lon"]) for r in rows if r["status"].startswith("assigned_")}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(13.5, 10), facecolor="#f8faf9")
    grid = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], width_ratios=[1.1, 1], left=.13, right=.96, top=.82, bottom=.155, hspace=.32, wspace=.32)
    ax = fig.add_subplot(grid[0, :])
    ax.set_facecolor("#e7eff1")
    land = json.loads(Path(args.land_geojson).read_text())
    for feature in land["features"]:
        geom = feature["geometry"]
        polygons = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
        for poly in polygons:
            coords = np.asarray(poly[0])
            ax.add_patch(Polygon(coords, facecolor="#d9dfd9", edgecolor="#b5bfb8", linewidth=.25))
    for category, marker in [("missing", "x"), ("coastal", "o"), ("cardinal", "o"), ("native", "o")]:
        pts = [point for point, states in sites.items() if category in states]
        if pts:
            ax.scatter([p[1] for p in pts], [p[0] for p in pts], c=COLORS[category], marker=marker, s=23 if marker == "x" else 32, linewidths=1, zorder=4)
    ax.set(xlim=(-112, 126), ylim=(-21, 37), xlabel="Longitud", ylabel="Latitud")
    ax.grid(alpha=.2, linewidth=.5)
    ax.set_title("Distribución de las coordenadas disponibles", loc="left", fontsize=12, fontweight="bold", pad=10)
    labels = {"native": "Celda nativa", "coastal": "Revisión costera previa", "cardinal": "Búsqueda cardinal", "missing": "Sin valor asignado"}
    ax.legend(handles=[Patch(facecolor=COLORS[k], label=labels[k]) for k in COLORS if counts[k]], loc="lower right", ncol=2, framealpha=.96, fontsize=9)
    bar = fig.add_subplot(grid[1, 0])
    ordered = sorted(countries, key=lambda c: (-sum(countries[c].values()), c))
    y = np.arange(len(ordered))
    left = np.zeros(len(ordered))
    for category in COLORS:
        values = np.array([countries[c][category] for c in ordered])
        bar.barh(y, values, left=left, color=COLORS[category], height=.66)
        left += values
    for yi, country, total in zip(y, ordered, left):
        bar.text(total+2, yi, f"{int(total)-countries[country]['missing']}/{int(total)}", va="center", fontsize=9)
    bar.set(yticks=y, yticklabels=ordered, xlabel="Registros; etiqueta = con valor / total", xlim=(0, max(left)*1.2))
    bar.invert_yaxis()
    bar.spines['left'].set_visible(False)
    bar.set_title("Cobertura por país o territorio", loc="left", fontsize=12, fontweight="bold", pad=12)
    note = fig.add_subplot(grid[1, 1]);note.axis("off")
    note.text(0, 1, "Método de recuperación", fontsize=14, fontweight="bold", va="top", color="#163e46")
    status = collections.Counter(r["status"] for r in rows)
    lines = [f"{counts['native']}  registros en celdas nativas", f"{counts['coastal']}  en celdas costeras revisadas previamente"]
    if counts['cardinal']:
        lines.append(f"{counts['cardinal']}  recuperados por búsqueda N / E / S / O")
    lines += ["", f"{counts['missing']}  registros conservados con NA:"]
    descriptions = {"held_by_spatial_review": "retenidos tras revisión espacial", "held_existing_coordinate_warning": "con advertencias previas de coordenadas", "no_native_coverage": "sin cobertura nativa ni alternativa aprobada", "no_cardinal_coverage_within_limit": "sin dato en la búsqueda cardinal hasta el límite", "missing_or_invalid_coordinates": "sin coordenadas válidas"}
    lines += [f"  {status[s]}  {label}" for s, label in descriptions.items() if status[s]]
    lines += ["", f"{len(selected_sites)}  coordenadas con valores asignados", f"{len(cells)}  celdas marinas distintas utilizadas"]
    if 'previous_status' in rows[0]:
        lines += ["", "Las advertencias previas siguen en la auditoría."]
    text = "\n".join(lines)
    note.text(0, .86, text, fontsize=11, va="top", linespacing=1.45, color="#263c42")
    fig.text(.07, .956, "Salinidad marina de referencia", fontsize=25, fontweight="bold", color="#123a44")
    fig.text(.07, .918, "Bio-ORACLE v3  ·  2000–2019  ·  Media y amplitud anual media  ·  Resolución 0.05°", fontsize=12, color="#52666b")
    fig.text(.07, .868, f"{len(rows)} registros conservados     |     {assigned} con valores     |     {counts['missing']} con NA", fontsize=17, fontweight="bold", color="#187b83")
    fig.text(.07, .048, "Contexto marino superficial; no representa salinidad intersticial ni la fecha de colecta.\nUn punto cartográfico puede representar varios registros. Fuente cartográfica: Natural Earth (dominio público).", fontsize=10, color="#52666b", linespacing=1.5)
    prefix = Path(args.out_prefix);prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(prefix)+".png", dpi=180)
    fig.savefig(str(prefix)+".svg")
    plt.close(fig)
    provenance = {"audit_sha256": hashlib.sha256(Path(args.audit).read_bytes()).hexdigest(), "land_sha256": hashlib.sha256(Path(args.land_geojson).read_bytes()).hexdigest(), "n_rows": len(rows), "n_assigned": assigned, "n_mapped_coordinate_pairs": len(sites), "n_assigned_coordinate_pairs": len(selected_sites), "n_assigned_cells": len(cells)}
    Path(str(prefix)+".json").write_text(json.dumps(provenance, indent=2)+"\n")
    print(json.dumps(provenance))


if __name__ == "__main__":
    main()
