#!/usr/bin/env python3
"""Make leakage-safe spatial folds for the development graph."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path

def read(p):
    with p.open(newline="") as f: return list(csv.DictReader(f))

def run(root: Path, out: Path, folds: int = 3) -> None:
    nodes = read(root / "slice001_grain_nodes_clean.csv")
    edges = read(root / "slice001_grain_boundaries_clean.csv")
    # 2D development split; 3D production split must group all voxels of a physical grain.
    xs = sorted(float(n["centroid_x_um"]) for n in nodes)
    cuts = [xs[min(len(xs) - 1, (i + 1) * len(xs) // folds)] for i in range(folds - 1)]
    def fold(x): return sum(x >= c for c in cuts)
    node_fold = {int(n["grain_id"]): fold(float(n["centroid_x_um"])) for n in nodes}
    result = {"folds": []}
    for k in range(folds):
        test_nodes = sorted(g for g, f in node_fold.items() if f == k)
        train_nodes = sorted(g for g, f in node_fold.items() if f != k)
        train_set = set(train_nodes)
        train_edges = [[int(e["grain_a"]), int(e["grain_b"])] for e in edges
                       if int(e["grain_a"]) in train_set and int(e["grain_b"]) in train_set]
        result["folds"].append({"fold": k, "train_nodes": train_nodes,
                                "test_nodes": test_nodes, "train_edges": train_edges})
    out.mkdir(parents=True, exist_ok=True)
    (out / "spatial_folds.json").write_text(json.dumps(result, indent=2))
    print(f"folds={folds} output={out / 'spatial_folds.json'}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--root", type=Path, required=True); p.add_argument("--out", type=Path, required=True); p.add_argument("--folds", type=int, default=3); a = p.parse_args(); run(a.root, a.out, a.folds)
