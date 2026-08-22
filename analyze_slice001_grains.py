#!/usr/bin/env python3
"""Approximate grain and grain-boundary chemistry from one EBSD slice."""
from __future__ import annotations

import csv
import math
import itertools
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
NC, NR = 549, 478
N = NC * NR
MISORIENTATION_DEG = float(os.getenv("A617_MISORIENTATION_DEG", "5"))
ELEMENTS = ("o", "cr", "ni", "fe")


class DSU:
    def __init__(self, n: int) -> None:
        self.p = list(range(n))

    def find(self, x: int) -> int:
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a: int, b: int) -> None:
        a, b = self.find(a), self.find(b)
        if a != b:
            self.p[b] = a


def angle_diff(a: float, b: float) -> float:
    d = abs(a - b)
    return min(d, 2 * math.pi - d)


def mat_quat(m: tuple[tuple[int, int, int], ...]) -> tuple[float, float, float, float]:
    t = 1 + m[0][0] + m[1][1] + m[2][2]
    if t > 0:
        s = 2 * math.sqrt(t)
        return ((s / 4), (m[2][1] - m[1][2]) / s, (m[0][2] - m[2][0]) / s, (m[1][0] - m[0][1]) / s)
    i = max(range(3), key=lambda k: m[k][k])
    j, k = (i + 1) % 3, (i + 2) % 3
    s = 2 * math.sqrt(max(1e-12, 1 + m[i][i] - m[j][j] - m[k][k]))
    q = [0.0, 0.0, 0.0, 0.0]
    q[i + 1] = s / 4
    q[0] = (m[k][j] - m[j][k]) / s
    q[j + 1] = (m[j][i] + m[i][j]) / s
    q[k + 1] = (m[k][i] + m[i][k]) / s
    return tuple(q)


def cubic_symmetry() -> list[tuple[float, float, float, float]]:
    out = []
    for perm in itertools.permutations(range(3)):
        inv = sum(perm[i] > perm[j] for i in range(3) for j in range(i + 1, 3))
        for signs in itertools.product((-1, 1), repeat=3):
            if ((-1) ** inv) * signs[0] * signs[1] * signs[2] == 1:
                m = tuple(tuple(signs[i] if perm[i] == j else 0 for j in range(3)) for i in range(3))
                out.append(mat_quat(m))
    assert len(out) == 24
    return out


def euler_quat(a: float, b: float, c: float) -> tuple[float, float, float, float]:
    ca, sa, cb, sb, cc, sc = math.cos(a), math.sin(a), math.cos(b), math.sin(b), math.cos(c), math.sin(c)
    m = ((ca * cc - sa * sc * cb, sa * cc + ca * sc * cb, sc * sb),
         (-ca * sc - sa * cc * cb, -sa * sc + ca * cc * cb, cc * sb),
         (sa * sb, -ca * sb, cb))
    return mat_quat(m)


def qmul(a, b):
    return (a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
            a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
            a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
            a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0])


def misorientation(q1, q2, syms) -> float:
    inv = (q1[0], -q1[1], -q1[2], -q1[3])
    best = 0.0
    for s in syms:
        q = qmul(s, qmul(q2, inv))
        best = max(best, abs(q[0]))
    return math.degrees(2 * math.acos(min(1.0, best)))


def load_points() -> tuple[list[tuple[float, float, float, int, float]], list[list[int | None]]]:
    ori: list[tuple[float, float, float, int, float]] = []
    with (ROOT / "slice001_ebsd_points.csv").open() as f:
        for row in csv.DictReader(f):
            ori.append((float(row["phi1"]), float(row["phi"]), float(row["phi2"]), int(float(row["phase_index"])), float(row["ci"])))
    assert len(ori) == N
    eds = [[None] * N for _ in ELEMENTS]
    with (ROOT / "slice001_eds_features.csv").open() as f:
        for i, row in enumerate(csv.DictReader(f)):
            for j, e in enumerate(ELEMENTS):
                v = row[f"eds_{e}k_counts"]
                eds[j][i] = None if not v else int(v)
    return ori, eds


def main() -> None:
    ori, eds = load_points()
    syms = cubic_symmetry()
    quats = [euler_quat(a, b, c) for a, b, c, _, _ in ori]
    dsu = DSU(N)
    for r in range(NR):
        for c in range(NC):
            i = r * NC + c
            for j in ((i + 1) if c + 1 < NC else -1, (i + NC) if r + 1 < NR else -1):
                if j < 0 or ori[i][3] != ori[j][3]:
                    continue
                d = misorientation(quats[i], quats[j], syms)
                if d <= MISORIENTATION_DEG:
                    dsu.union(i, j)

    roots = {}
    labels = []
    for i in range(N):
        root = dsu.find(i)
        labels.append(roots.setdefault(root, len(roots)))
    grain_count = len(roots)

    boundary_pixels = [False] * N
    pairs: dict[tuple[int, int], list[int]] = {}
    for r in range(NR):
        for c in range(NC):
            i = r * NC + c
            for j in ((i + 1) if c + 1 < NC else -1, (i + NC) if r + 1 < NR else -1):
                if j < 0 or labels[i] == labels[j]:
                    continue
                boundary_pixels[i] = boundary_pixels[j] = True
                pair = tuple(sorted((labels[i], labels[j])))
                pairs.setdefault(pair, []).extend((i, j))

    sizes = [0] * grain_count
    sums = [[0.0, 0.0] for _ in range(grain_count)]
    interior_sum = [[0.0] * len(ELEMENTS) for _ in range(grain_count)]
    interior_count = [[0] * len(ELEMENTS) for _ in range(grain_count)]
    for i, g in enumerate(labels):
        sizes[g] += 1
        sums[g][0] += ori[i][0] * 0 + (i % NC) * 0.1
        sums[g][1] += (i // NC) * 0.0866
        if not boundary_pixels[i]:
            for k in range(len(ELEMENTS)):
                if eds[k][i] is not None:
                    interior_sum[g][k] += eds[k][i]
                    interior_count[g][k] += 1

    with (ROOT / "slice001_grain_nodes.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["grain_id", "pixel_count", "area_um2", "centroid_x_um", "centroid_y_um"])
        for g in range(grain_count):
            w.writerow([g, sizes[g], sizes[g] * 0.1 * 0.0866, sums[g][0] / sizes[g], sums[g][1] / sizes[g]])

    with (ROOT / "slice001_grain_labels.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["row", "col", "grain_id", "boundary_pixel"])
        for i, g in enumerate(labels):
            w.writerow([i // NC, i % NC, g, int(boundary_pixels[i])])

    with (ROOT / "slice001_grain_boundaries.csv").open("w", newline="") as f:
        w = csv.writer(f)
        header = ["grain_a", "grain_b", "boundary_pixel_count"]
        for e in ELEMENTS:
            header += [f"{e}_boundary_mean", f"{e}_bulk_mean", f"{e}_relative_enrichment", f"{e}_depletion"]
        w.writerow(header)
        for (a, b), ids in sorted(pairs.items()):
            row = [a, b, len(ids) // 2]
            for k in range(len(ELEMENTS)):
                boundary = [eds[k][i] for i in ids if eds[k][i] is not None]
                bm = sum(boundary) / len(boundary) if boundary else float("nan")
                n_bulk = interior_count[a][k] + interior_count[b][k]
                im = (interior_sum[a][k] + interior_sum[b][k]) / n_bulk if n_bulk else float("nan")
                rel = bm / im if im and not math.isnan(bm) else float("nan")
                dep = (im - bm) / im if im and not math.isnan(bm) else float("nan")
                row += [bm, im, rel, dep]
            w.writerow(row)
    print(f"grains={grain_count} boundaries={len(pairs)} threshold_deg={MISORIENTATION_DEG}")


if __name__ == "__main__":
    main()
