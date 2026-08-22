#!/usr/bin/env python3
"""Build a small, reproducible index for the downloaded Slice 001 data."""
from __future__ import annotations

import csv
import re
import struct
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ANG = ROOT / "sqr_ebsd-sliceimage-001_mod.ang"


def png_size(path: Path) -> tuple[int, int] | None:
    with path.open("rb") as f:
        if f.read(8) != b"\x89PNG\r\n\x1a\n":
            return None
        length = struct.unpack(">I", f.read(4))[0]
        if f.read(4) != b"IHDR" or length < 8:
            return None
        width, height = struct.unpack(">II", f.read(8))
        return width, height


def parse_ang(path: Path) -> tuple[dict[str, str], list[str], list[list[str]]]:
    meta: dict[str, str] = {}
    columns: list[str] = []
    rows: list[list[str]] = []
    data_started = False
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            if line.startswith("# COLUMN_HEADERS:"):
                columns = [x.strip().lower().replace(" ", "_") for x in line.split(":", 1)[1].split(",")]
            else:
                m = re.match(r"#\s*([^:]+):\s*(.*)", line)
                if m:
                    meta[m.group(1).strip().lower()] = m.group(2).strip()
            continue
        data_started = True
        if data_started:
            values = line.split()
            if columns and len(values) >= 10:
                rows.append(values)
    if not columns or not rows:
        raise RuntimeError(f"Could not parse ANG data: {path}")
    return meta, columns, rows


def classify(name: str) -> str:
    lower = name.lower()
    if lower.endswith(".ang"):
        return "ebsd_text"
    if lower.endswith(".spd") or lower.endswith(".dat") or lower.endswith(".ipr"):
        return "binary_or_instrument"
    if "sem-image" in lower:
        return "sem_image"
    if lower.endswith(".tif") or lower.endswith(".tiff"):
        return "element_image_tiff"
    if lower.endswith(".png"):
        return "element_or_ebsd_image"
    return "other"


def read_eds(path: Path, ncols: int, nrows: int) -> list[int | None]:
    raw = path.read_bytes()
    assert len(raw) == 36 + ncols * nrows * 4, (path.name, len(raw))
    header = struct.unpack("<9I", raw[:36])
    assert header[:2] == (ncols, nrows), (path.name, header[:2])
    # The vendor file stores five footer words after the pixel grid.
    values = list(struct.unpack(f"<{ncols * nrows}I", raw[36:]))
    return [None if x == 0x80000000 else x for x in values[:-5]] + [None] * 5


def main() -> None:
    meta, columns, rows = parse_ang(ANG)
    ncols = int(meta["ncols_odd"])
    nrows = int(meta["nrows"])
    assert len(rows) == ncols * nrows, (len(rows), ncols, nrows)
    with (ROOT / "slice001_ebsd_points.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["row", "col", "x_um", "y_um", *columns])
        for i, values in enumerate(rows):
            padded = values + [""] * max(0, len(columns) - len(values))
            writer.writerow([i // ncols, i % ncols, values[3], values[4], *padded[: len(columns)]])

    eds_files = sorted(ROOT.glob("ebsd-sliceimage-001_0_*.dat"))
    eds_files = [p for p in eds_files if not p.name.endswith("_cps.dat")]
    eds = {p.stem.rsplit("_", 1)[-1].upper(): read_eds(p, ncols, nrows) for p in eds_files}
    with (ROOT / "slice001_eds_features.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        channels = sorted(eds)
        writer.writerow(["row", "col", *[f"eds_{x.lower()}_counts" for x in channels]])
        for i in range(ncols * nrows):
            writer.writerow([i // ncols, i % ncols, *("" if eds[x][i] is None else eds[x][i] for x in channels)])

    with (ROOT / "slice001_manifest.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "bytes", "category", "width_px", "height_px", "orientation_transform", "coordinate_source"])
        for path in sorted(ROOT.iterdir()):
            if not path.is_file() or path.name in {"slice001_manifest.csv", "slice001_ebsd_points.csv"}:
                continue
            size = png_size(path) if path.suffix.lower() == ".png" else None
            category = classify(path.name)
            # Slice 001 follows the dominant EBSD direction; preserve source files and record
            # this canonical-orientation assumption for later outlier review.
            transform = "rotate_180" if path.name == "ebsd-sliceimage-001_ipf-iq.png" else "identity"
            source = "manual_flip_vertical_horizontal" if transform == "rotate_180" else ("image_export_majority_direction" if category in {"element_or_ebsd_image", "element_image_tiff", "sem_image"} else "ebsd_grid")
            writer.writerow([path.name, path.stat().st_size, category, *(size or ("", "")), transform, source])

    print(f"parsed_rows={len(rows)} grid={ncols}x{nrows}")
    print(f"columns={len(columns)} eds_channels={len(eds)} manifest=slice001_manifest.csv points=slice001_ebsd_points.csv eds=slice001_eds_features.csv")


if __name__ == "__main__":
    main()
