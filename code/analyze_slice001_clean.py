#!/usr/bin/env python3
"""EBSD cleanup, grain reconstruction, and boundary chemistry for Slice 001."""
from __future__ import annotations
import csv, math, os
from pathlib import Path
from analyze_slice001_grains import (NC, NR, N, ELEMENTS, DSU, cubic_symmetry,
    euler_quat, misorientation, load_points)

ROOT = Path(__file__).resolve().parent
CI_MIN = float(os.getenv("A617_CI_MIN", "0.1"))
MIN_PIXELS = int(os.getenv("A617_MIN_GRAIN_PIXELS", "20"))
THRESHOLD = float(os.getenv("A617_MISORIENTATION_DEG", "5"))

def neighbours(i):
    r, c = divmod(i, NC)
    return ([i-NC] if r else []) + ([i+NC] if r+1 < NR else []) + ([i-1] if c else []) + ([i+1] if c+1 < NC else [])

def fill_invalid(labels):
    for _ in range(4):
        for i, g in enumerate(labels):
            if g >= 0: continue
            ns = [labels[j] for j in neighbours(i) if labels[j] >= 0]
            if ns: labels[i] = max(set(ns), key=ns.count)

def main():
    ori, eds = load_points()
    valid = [x[4] >= CI_MIN and x[3] > 0 for x in ori]
    quats = [euler_quat(x[0], x[1], x[2]) for x in ori]
    syms = cubic_symmetry(); dsu = DSU(N)
    for i in range(N):
        if not valid[i]: continue
        for j in (i + 1 if i % NC < NC - 1 else -1, i + NC if i // NC < NR - 1 else -1):
            if j >= 0 and valid[j] and ori[i][3] == ori[j][3] and misorientation(quats[i], quats[j], syms) <= THRESHOLD:
                dsu.union(i, j)
    roots, labels = {}, []
    for i in range(N):
        if not valid[i]: labels.append(-1); continue
        root = dsu.find(i); labels.append(roots.setdefault(root, len(roots)))
    sizes = [0] * len(roots)
    for g in labels:
        if g >= 0: sizes[g] += 1
    labels = [g if g >= 0 and sizes[g] >= MIN_PIXELS else -1 for g in labels]
    fill_invalid(labels)
    remap = {}
    labels = [remap.setdefault(g, len(remap)) if g >= 0 else -1 for g in labels]
    ng = len(remap)

    boundary = [False] * N; pairs = {}
    for i in range(N):
        if labels[i] < 0: continue
        for j in (i + 1 if i % NC < NC - 1 else -1, i + NC if i // NC < NR - 1 else -1):
            if j < 0 or labels[j] < 0 or labels[i] == labels[j]: continue
            boundary[i] = boundary[j] = True
            pairs.setdefault(tuple(sorted((labels[i], labels[j]))), []).extend((i, j))

    size = [0] * ng; sx = [0.0] * ng; sy = [0.0] * ng
    isum = [[0.0] * 4 for _ in range(ng)]; icount = [[0] * 4 for _ in range(ng)]
    for i, g in enumerate(labels):
        if g < 0: continue
        size[g] += 1; sx[g] += (i % NC) * 0.1; sy[g] += (i // NC) * 0.0866
        if not boundary[i]:
            for k in range(4):
                if eds[k][i] is not None: isum[g][k] += eds[k][i]; icount[g][k] += 1

    with (ROOT / "slice001_grain_nodes_clean.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["grain_id", "pixel_count", "area_um2", "centroid_x_um", "centroid_y_um"])
        for g in range(ng): w.writerow([g, size[g], size[g]*0.1*0.0866, sx[g]/size[g], sy[g]/size[g]])
    with (ROOT / "slice001_grain_labels_clean.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["row", "col", "grain_id", "boundary_pixel"])
        for i, g in enumerate(labels): w.writerow([i//NC, i%NC, g, int(boundary[i])])
    with (ROOT / "slice001_grain_boundaries_clean.csv").open("w", newline="") as f:
        h = ["grain_a", "grain_b", "boundary_pixel_count"]
        for e in ELEMENTS: h += [f"{e}_boundary_mean", f"{e}_bulk_mean", f"{e}_relative_enrichment", f"{e}_depletion"]
        w = csv.writer(f); w.writerow(h)
        for (a, b), ids in sorted(pairs.items()):
            row = [a, b, len(ids)//2]
            for k in range(4):
                vals = [eds[k][i] for i in ids if eds[k][i] is not None]
                bm = sum(vals)/len(vals) if vals else float("nan")
                n = icount[a][k] + icount[b][k]; im = (isum[a][k]+isum[b][k])/n if n else float("nan")
                row += [bm, im, bm/im if im and not math.isnan(bm) else float("nan"), (im-bm)/im if im and not math.isnan(bm) else float("nan")]
            w.writerow(row)
    print(f"valid={sum(valid)} grains={ng} boundaries={len(pairs)} ci_min={CI_MIN} min_pixels={MIN_PIXELS} threshold_deg={THRESHOLD}")

if __name__ == "__main__": main()
