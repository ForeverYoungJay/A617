#!/usr/bin/env python3
"""Build the leakage-controlled Slice 001 edge-learning contract."""
from __future__ import annotations
import argparse, csv, math
from pathlib import Path

ELEMENTS = ("o", "cr", "ni", "fe")

def read(path):
    with path.open(newline="") as f: return list(csv.DictReader(f))

def f(row, key, default=0.0):
    try: return float(row[key])
    except (KeyError, TypeError, ValueError): return default

def logratio(boundary, bulk):
    eps = 1e-9
    return math.log((boundary + eps) / (bulk + eps)) if boundary >= 0 and bulk >= 0 else math.nan

def run(root: Path, out: Path, include_reconstructed: bool = False) -> None:
    nodes = read(root / "slice001_grain_nodes_clean.csv")
    edges = read(root / "slice001_grain_boundaries_clean.csv")
    labels = read(root / "slice001_grain_labels_clean.csv")
    # Node features intentionally exclude IQ/CI, gap, absolute centroid and chemistry.
    degree = {int(n["grain_id"]): 0 for n in nodes}
    for e in edges:
        a, b = int(e["grain_a"]), int(e["grain_b"]); degree[a] += 1; degree[b] += 1
    node_out = out / "phase1_node_features.csv"
    with node_out.open("w", newline="") as g:
        fields = ["grain_id", "area_um2", "pixel_count", "degree"]
        w = csv.DictWriter(g, fieldnames=fields); w.writeheader()
        for n in nodes:
            w.writerow({"grain_id": n["grain_id"], "area_um2": n["area_um2"],
                        "pixel_count": n["pixel_count"], "degree": degree[int(n["grain_id"])]})

    # Existing clean boundaries are direct contacts. Reconstructed edges live in the
    # separate confidence table and are deliberately excluded from the primary set.
    reconstructed = {(int(r["grain_a"]), int(r["grain_b"]))
                     for r in read(root / "slice001_reconstructed_edges_confidence.csv")}
    fields = ["grain_a", "grain_b", "boundary_pixel_count", "edge_set", "label_quality"]
    for e in ELEMENTS: fields.append(f"logratio_{e}")
    edge_out = out / ("phase1_edges_all.csv" if include_reconstructed else "phase1_edges_direct.csv")
    with edge_out.open("w", newline="") as g:
        w = csv.DictWriter(g, fieldnames=fields); w.writeheader()
        for e in edges:
            pair = (int(e["grain_a"]), int(e["grain_b"]))
            if not include_reconstructed and pair in reconstructed: continue
            target = {f"logratio_{x}": logratio(f(e, f"{x}_boundary_mean"), f(e, f"{x}_bulk_mean")) for x in ELEMENTS}
            valid = [v for v in target.values() if not math.isnan(v)]
            w.writerow({"grain_a": pair[0], "grain_b": pair[1],
                        "boundary_pixel_count": e["boundary_pixel_count"],
                        "edge_set": "reconstructed" if pair in reconstructed else "direct",
                        "label_quality": len(valid) / len(ELEMENTS), **target})
    print(f"nodes={len(nodes)} edges={len(edges)} output={edge_out}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--root", type=Path, required=True); p.add_argument("--out", type=Path, required=True)
    p.add_argument("--include-reconstructed", action="store_true"); a = p.parse_args(); a.out.mkdir(parents=True, exist_ok=True); run(a.root, a.out, a.include_reconstructed)
