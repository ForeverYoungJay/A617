#!/usr/bin/env python3
"""Validate reconstructed GB geometry, topology, and parameter stability."""
from __future__ import annotations
import csv, math
from collections import Counter
from pathlib import Path
from analyze_slice001_clean import NC, NR, N, neighbours, fill_invalid
from analyze_slice001_grains import DSU, cubic_symmetry, euler_quat, load_points, misorientation

ROOT = Path(__file__).resolve().parent
DX, DY = 0.1, 0.0866

def build(threshold: float, min_pixels: int):
    ori, _ = load_points(); valid = [x[4] >= 0.1 and x[3] > 0 for x in ori]
    q = [euler_quat(x[0], x[1], x[2]) for x in ori]; syms = cubic_symmetry(); d = DSU(N)
    for i in range(N):
        if not valid[i]: continue
        for j in (i + 1 if i % NC < NC - 1 else -1, i + NC if i // NC < NR - 1 else -1):
            if j >= 0 and valid[j] and ori[i][3] == ori[j][3] and misorientation(q[i], q[j], syms) <= threshold: d.union(i, j)
    roots, labels = {}, []
    for i in range(N):
        if not valid[i]: labels.append(-1); continue
        r = d.find(i); labels.append(roots.setdefault(r, len(roots)))
    sizes = Counter(x for x in labels if x >= 0)
    labels = [x if x >= 0 and sizes[x] >= min_pixels else -1 for x in labels]
    raw = labels[:]; fill_invalid(labels)
    remap = {}; labels = [remap.setdefault(x, len(remap)) if x >= 0 else -1 for x in labels]
    raw_remap = {old: new for new, old in enumerate(sorted(set(x for x in raw if x >= 0)))}
    raw = [raw_remap.get(x, -1) for x in raw]
    return raw, labels, len(remap)

def edges(labels, raw):
    out = {}
    for i in range(N):
        if labels[i] < 0: continue
        for j, length in ((i + 1 if i % NC < NC - 1 else -1, DY), (i + NC if i // NC < NR - 1 else -1, DX)):
            if j < 0 or labels[j] < 0 or labels[i] == labels[j]: continue
            p = tuple(sorted((labels[i], labels[j]))); x = out.setdefault(p, [0, 0.0, 0])
            x[0] += 1; x[1] += length
            if raw[i] >= 0 and raw[j] >= 0 and raw[i] != raw[j]: x[2] += 1
    return out

def main():
    raw, labels, ng = build(5, 20); es = edges(labels, raw)
    with (ROOT / "slice001_gb_validation.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["grain_a", "grain_b", "shared_pixel_contacts", "shared_length_um", "direct_valid_contacts", "edge_class"])
        for (a, b), (n, length, direct) in sorted(es.items()): w.writerow([a, b, n, length, direct, "direct" if direct else "reconstructed_across_mask"])
    deg = Counter();
    for a, b in es: deg[a] += 1; deg[b] += 1
    with (ROOT / "slice001_gb_sensitivity.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["threshold_deg", "min_pixels", "grain_count", "edge_count", "direct_edge_count", "mean_degree", "planar_bound"])
        for t in (3, 5, 7, 10):
            for m in (10, 20, 50, 100):
                r, l, n = build(t, m); e = edges(l, r); direct = sum(any(x[2] for x in [v]) for v in e.values()); w.writerow([t, m, n, len(e), direct, 2*len(e)/n if n else 0, 3*n-6])
    print(f"grains={ng} edges={len(es)} direct_edges={sum(v[2] > 0 for v in es.values())} reconstructed={sum(v[2] == 0 for v in es.values())} mean_degree={2*len(es)/ng:.3f}")

if __name__ == "__main__": main()
