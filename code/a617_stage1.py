#!/usr/bin/env python3
"""Server-side Slice 001+ ingestion for the Alloy 617 graph-data pipeline.

This stage deliberately stops at a lossless, coordinate-checked EBSD/EDS table.
Grain reconstruction is a separate stage because it needs the validated symmetry
and segmentation parameters, not a silent fallback heuristic.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import struct
from pathlib import Path


SENTINEL = 0x80000000


def ang(path: Path):
    meta, columns, rows = {}, [], []
    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            if line.startswith("# COLUMN_HEADERS:"):
                columns = [x.strip().lower().replace(" ", "_")
                           for x in line.split(":", 1)[1].split(",")]
            else:
                m = re.match(r"#\s*([^:]+):\s*(.*)", line)
                if m:
                    meta[m.group(1).strip().lower()] = m.group(2).strip()
            continue
        values = line.split()
        if len(values) >= 10:
            rows.append(values)
    ncols = int(meta.get("ncols_odd", meta.get("ncols", "0")))
    nrows = int(meta.get("nrows", "0"))
    if not columns or not ncols or not nrows or len(rows) != ncols * nrows:
        raise ValueError(f"invalid ANG grid: {path} ({len(rows)} rows, {ncols}x{nrows})")
    return ncols, nrows, columns, rows


def eds(path: Path, ncols: int, nrows: int) -> list[int | None]:
    raw = path.read_bytes()
    expected = 36 + ncols * nrows * 4
    if len(raw) not in (expected, expected + 20):
        raise ValueError(f"unexpected EDS size: {path.name}: {len(raw)}")
    h = struct.unpack("<9I", raw[:36])
    if h[:2] != (ncols, nrows):
        raise ValueError(f"EDS/ANG grid mismatch: {path.name}: {h[:2]} vs {(ncols, nrows)}")
    values = struct.unpack(f"<{ncols*nrows}I", raw[36:expected])
    return [None if x == SENTINEL else x for x in values]


def median(values):
    values = sorted(values)
    if not values:
        return math.nan
    n = len(values)
    return values[n // 2] if n % 2 else (values[n // 2 - 1] + values[n // 2]) / 2


def run(root: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    ang_files = sorted(root.glob("*.ang"))
    if not ang_files:
        raise SystemExit(f"no .ang files under {root}")
    manifest = []
    for ap in ang_files:
        ncols, nrows, columns, rows = ang(ap)
        stem = ap.stem.replace("sqr_", "").replace("_mod", "")
        eds_prefix = stem if "_0" in stem else f"{stem}_0"
        candidates = sorted(root.glob(f"{eds_prefix}_*.dat"))
        channels = {}
        for dp in candidates:
            name = dp.stem.rsplit("_", 1)[-1].lower()
            if name.endswith("_cps") or name in {"ebsd", "image"}:
                continue
            try:
                channels[name] = eds(dp, ncols, nrows)
            except ValueError:
                continue
        table = out / f"{stem}_features.csv"
        with table.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["slice", "row", "col", *columns,
                        *[f"eds_{k}_counts" for k in sorted(channels)]])
            for i, row in enumerate(rows):
                w.writerow([stem, i // ncols, i % ncols, *row,
                            *[("" if channels[k][i] is None else channels[k][i])
                              for k in sorted(channels)]])
        summary = {k: median([x for x in v if x is not None]) for k, v in channels.items()}
        manifest.append([stem, ap.name, ncols, nrows, len(rows), len(channels), table.name,
                         ";".join(f"{k}:{summary[k]:.6g}" for k in sorted(summary))])
    with (out / "slice_manifest.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["slice", "ang_file", "ncols", "nrows", "points", "eds_channels",
                    "feature_file", "channel_medians"])
        w.writerows(manifest)
    print(f"slices={len(manifest)} output={out}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True, help="server directory containing ANG and EDS DAT files")
    p.add_argument("--out", type=Path, required=True, help="output directory")
    args = p.parse_args()
    run(args.root, args.out)
