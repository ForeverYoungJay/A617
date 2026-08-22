#!/usr/bin/env python3
"""Check whether GB enrichment survives wider spatial sampling windows."""
from __future__ import annotations
import csv, math, statistics
from collections import defaultdict
from pathlib import Path
from validate_gb_slice001 import build
from analyze_slice001_grains import NC, NR, load_points

ROOT = Path(__file__).resolve().parent
ELEMENTS = ("o", "cr", "ni", "fe")
RADII = (1, 2, 3, 5)

def main():
    raw, labels, _ = build(5, 20)
    _, eds = load_points()
    contacts = defaultdict(set)
    for i in range(NC * NR):
        if labels[i] < 0: continue
        for j in (i + 1 if i % NC < NC - 1 else -1, i + NC if i // NC < NR - 1 else -1):
            if j >= 0 and labels[j] >= 0 and labels[i] != labels[j]:
                contacts[tuple(sorted((labels[i], labels[j])))].update((i, j))
    out = []
    for (a, b), edge in sorted(contacts.items()):
        for radius in RADII:
            near = set()
            for i in edge:
                r, c = divmod(i, NC)
                for rr in range(max(0, r-radius), min(NR, r+radius+1)):
                    for cc in range(max(0, c-radius), min(NC, c+radius+1)):
                        j = rr*NC + cc
                        if labels[j] in (a, b): near.add(j)
            bulk = {i for i in range(NC*NR) if labels[i] in (a, b) and i not in near}
            for k, e in enumerate(ELEMENTS):
                bv = [eds[k][i] for i in near if eds[k][i] is not None]
                iv = [eds[k][i] for i in bulk if eds[k][i] is not None]
                bm = sum(bv)/len(bv) if bv else float("nan"); im = sum(iv)/len(iv) if iv else float("nan")
                out.append([a,b,radius,e,len(bv),len(iv),bm,im,bm/im if im else float("nan")])
    with (ROOT / "slice001_eds_scale_robustness.csv").open("w", newline="") as f:
        w=csv.writer(f); w.writerow(["grain_a","grain_b","radius_px","element","near_n","bulk_n","near_mean","bulk_mean","relative_enrichment"]); w.writerows(out)
    with (ROOT / "slice001_eds_scale_summary.csv").open("w", newline="") as f:
        w=csv.writer(f); w.writerow(["radius_px","element","edge_n","median_relative_enrichment","q25","q75"])
        for radius in RADII:
            for e in ELEMENTS:
                vals=[float(x[-1]) for x in out if x[2]==radius and x[3]==e and math.isfinite(float(x[-1]))]
                vals.sort(); q=lambda p: vals[min(len(vals)-1,int(p*(len(vals)-1)))]
                w.writerow([radius,e,len(vals),statistics.median(vals),q(.25),q(.75)] if vals else [radius,e,0,"","",""])
    print(f"edges={len(contacts)} rows={len(out)} radii={RADII}")

if __name__ == "__main__": main()
