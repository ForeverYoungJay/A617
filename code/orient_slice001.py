#!/usr/bin/env python3
"""Choose PNG orientation by correlation with the matching EDS grid."""
from __future__ import annotations

import csv
import struct
import zlib
from pathlib import Path

from prepare_slice001 import read_eds

ROOT = Path(__file__).resolve().parent
N_COLS, N_ROWS = 549, 478


def png_gray(path: Path) -> tuple[int, int, list[float]]:
    b = path.read_bytes()
    assert b[:8] == b"\x89PNG\r\n\x1a\n"
    pos, width, height, color, depth, chunks = 8, 0, 0, 0, 0, bytearray()
    while pos < len(b):
        n = struct.unpack(">I", b[pos : pos + 4])[0]
        kind, payload = b[pos + 4 : pos + 8], b[pos + 8 : pos + 8 + n]
        pos += n + 12
        if kind == b"IHDR":
            width, height, depth, color = struct.unpack(">IIBB", payload[:10])
            assert depth == 8 and color in (0, 2, 4, 6), (path, depth, color)
        elif kind == b"IDAT":
            chunks.extend(payload)
        elif kind == b"IEND":
            break
    raw = zlib.decompress(chunks)
    channels = {0: 1, 2: 3, 4: 2, 6: 4}[color]
    stride = width * channels
    rows: list[bytes] = []
    prev = bytearray(stride)
    p = 0
    for _ in range(height):
        filt = raw[p]
        cur = bytearray(raw[p + 1 : p + 1 + stride])
        p += stride + 1
        for i in range(stride):
            left = cur[i - channels] if i >= channels else 0
            up = prev[i]
            ul = prev[i - channels] if i >= channels else 0
            if filt == 1:
                cur[i] = (cur[i] + left) & 255
            elif filt == 2:
                cur[i] = (cur[i] + up) & 255
            elif filt == 3:
                cur[i] = (cur[i] + ((left + up) // 2)) & 255
            elif filt == 4:
                q = left + up - ul
                pa, pb, pc = abs(q - left), abs(q - up), abs(q - ul)
                cur[i] = (cur[i] + (left if pa <= pb and pa <= pc else up if pb <= pc else ul)) & 255
            else:
                assert filt == 0, filt
        rows.append(bytes(cur))
        prev = cur
    # Nearest-neighbour sampling preserves the vendor raster's spatial ordering.
    out: list[float] = []
    for r in range(N_ROWS):
        sy = min(height - 1, int((r + 0.5) * height / N_ROWS))
        row = rows[sy]
        for c in range(N_COLS):
            sx = min(width - 1, int((c + 0.5) * width / N_COLS)) * channels
            out.append(sum(row[sx : sx + min(3, channels)]) / min(3, channels))
    return width, height, out


def corr(a: list[float], b: list[int | None]) -> float:
    pairs = [(x, y) for x, y in zip(a, b) if y is not None]
    ma = sum(x for x, _ in pairs) / len(pairs)
    mb = sum(y for _, y in pairs) / len(pairs)
    num = sum((x - ma) * (y - mb) for x, y in pairs)
    da = sum((x - ma) ** 2 for x, _ in pairs) ** 0.5
    db = sum((y - mb) ** 2 for _, y in pairs) ** 0.5
    return num / (da * db) if da and db else 0.0


def transformed(values: list[float], mode: str) -> list[float]:
    out = [0.0] * len(values)
    for r in range(N_ROWS):
        for c in range(N_COLS):
            rr, cc = (N_ROWS - 1 - r if mode in ("flip_vertical", "rotate_180") else r,
                      N_COLS - 1 - c if mode in ("flip_horizontal", "rotate_180") else c)
            out[r * N_COLS + c] = values[rr * N_COLS + cc]
    return out


def main() -> None:
    modes = ("identity", "flip_vertical", "flip_horizontal", "rotate_180")
    results: dict[str, tuple[str, float, float, float, float]] = {}
    for png in sorted(ROOT.glob("ebsd-sliceimage-001_*.png")):
        suffix = png.stem.rsplit("_", 1)[-1].lower()
        dat = ROOT / f"ebsd-sliceimage-001_0_{suffix}k.dat"
        if not dat.exists() or suffix in {"allelements", "ipf-iq", "ipf", "iq", "phase", "video"}:
            continue
        _, _, pixels = png_gray(png)
        eds = read_eds(dat, N_COLS, N_ROWS)
        scores = {mode: corr(transformed(pixels, mode), eds) for mode in modes}
        ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        results[png.name] = (ranked[0][0], ranked[0][1], scores["identity"], ranked[1][1], ranked[0][1] - ranked[1][1])

    manifest = ROOT / "slice001_manifest.csv"
    rows = list(csv.DictReader(manifest.open()))
    for row in rows:
        if row["file"] in results:
            best, score, identity, second, margin = results[row["file"]]
            row["orientation_transform"] = best
            row["coordinate_source"] = f"png_vs_eds_corr:{score:.4f};margin:{margin:.4f}"
    with manifest.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    for name, result in results.items():
        print(name, "best=%s score=%.4f identity=%.4f second=%.4f margin=%.4f" % result)


if __name__ == "__main__":
    main()
