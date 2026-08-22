#!/usr/bin/env python3
"""Write reconstructed-edge confidence data and SVG overlays."""
from __future__ import annotations
import csv, struct, zlib
from collections import defaultdict
from pathlib import Path
from validate_gb_slice001 import build
from analyze_slice001_grains import load_points, NC, NR

ROOT = Path(__file__).resolve().parent
W, H = 1536, 1157
DX, DY = W / NC, H / NR

def components(points):
    todo = set(points); out = []
    while todo:
        seed = todo.pop(); comp = [seed]; stack = [seed]
        while stack:
            i = stack.pop(); r, c = divmod(i, NC)
            for rr in range(max(0, r-1), min(NR, r+2)):
                for cc in range(max(0, c-1), min(NC, c+2)):
                    j = rr * NC + cc
                    if j in todo:
                        todo.remove(j); stack.append(j); comp.append(j)
        out.append(comp)
    return out

def png_gray(values):
    rows = []
    for r in range(NR):
        row = bytearray()
        for c in range(NC):
            row.append(max(0, min(255, int(values[r * NC + c] * 255))))
        rows.append(b"\0" + bytes(row))
    raw = b"".join(rows)
    def chunk(k, p):
        return struct.pack(">I", len(p)) + k + p + struct.pack(">I", zlib.crc32(k + p) & 0xffffffff)
    return (b"\x89PNG\r\n\x1a\n" +
            chunk(b"IHDR", struct.pack(">IIBBBBB", NC, NR, 8, 0, 0, 0, 0)) +
            chunk(b"IDAT", zlib.compress(raw, 6)) + chunk(b"IEND", b""))

def main():
    raw, labels, _ = build(5, 20)
    ori, _ = load_points()
    contacts = defaultdict(list)
    for i in range(NR * NC):
        if labels[i] < 0:
            continue
        for j in (i + 1 if i % NC < NC - 1 else -1,
                  i + NC if i // NC < NR - 1 else -1):
            if j < 0 or labels[j] < 0 or labels[i] == labels[j]:
                continue
            p = tuple(sorted((labels[i], labels[j])))
            direct = raw[i] >= 0 and raw[j] >= 0 and raw[i] != raw[j]
            contacts[p].append((i, j, direct))
    confidence = []
    for edge_id, (pair, pts) in enumerate(sorted(contacts.items()), 1):
        direct = sum(x[2] for x in pts)
        gap = len(pts) - direct
        cis = [ori[i][4] for i, j, _ in pts] + [ori[j][4] for i, j, _ in pts]
        low = sum(x < 0.1 for x in cis) / len(cis)
        support = min(1.0, len(pts) / 10.0)
        ci_support = sum(cis) / len(cis)
        # ponytail: heuristic confidence; calibrate when an annotated edge set exists.
        score = support * (0.5 * ci_support + 0.5 * (1 - low))
        endpoints = {x for p in pts for x in p[:2]}
        comps = components(endpoints)
        spans = [max(max(x % NC for x in comp)-min(x % NC for x in comp), max(x // NC for x in comp)-min(x // NC for x in comp)) for comp in comps]
        length = sum(0.0866 if j == i + 1 else 0.1 for i, j, _ in pts)
        if direct:
            continue
        score /= max(1, len(comps))
        confidence.append((edge_id, pair, pts, length, direct, gap, low, score, len(comps), max(spans, default=0), comps))
    with (ROOT / "slice001_reconstructed_edges_confidence.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["edge_id", "grain_a", "grain_b", "contact_pixels", "shared_length_um",
                    "direct_contacts", "gap_contacts", "low_ci_fraction", "heuristic_confidence",
                    "contact_components", "max_component_span_px", "disconnected_contact"])
        for eid, (a, b), pts, length, direct, gap, low, score, ncomp, span, comps in confidence:
            w.writerow([eid, a, b, len(pts), length, direct, gap, low, score, ncomp, span, int(ncomp > 1)])

    ci = [min(1.0, max(0.0, x[4])) for x in ori]
    (ROOT / "slice001_ci.png").write_bytes(png_gray(ci))
    review_ids = {12, 26, 35, 63, 83}
    for name in ("ebsd-sliceimage-001_ipf-iq.png", "ebsd-sliceimage-001_iq.png",
                 "ebsd-sliceimage-001_phase.png", "slice001_ci.png"):
        image_tag = (f'<image href="{name}" x="0" y="0" width="{W}" height="{H}" transform="translate({W} {H}) scale(-1 -1)"/>'
                     if name == "ebsd-sliceimage-001_ipf-iq.png" else
                     f'<image href="{name}" x="0" y="0" width="{W}" height="{H}"/>')
        svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">', image_tag]
        for eid, (a, b), pts, length, direct, gap, low, score, ncomp, span, comps in confidence:
            # Draw only local pixel-interface segments. Never connect distant contacts.
            for ci, comp in enumerate(comps):
                comp_set = set(comp); segments = []
                for i, j, _ in pts:
                    if i not in comp_set and j not in comp_set:
                        continue
                    x1, y1 = (i % NC + 0.5) * DX, (i // NC + 0.5) * DY
                    x2, y2 = (j % NC + 0.5) * DX, (j // NC + 0.5) * DY
                    if abs(x1 - x2) > DX * 1.1:
                        continue
                    segments.append((x1, y1, x2, y2))
                for x1, y1, x2, y2 in segments:
                    color = "#00ffff" if eid in review_ids else "#ff2020"
                    svg.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" stroke-width="2"/>')
                if comp and eid in review_ids:
                    cx = sum(i % NC for i in comp) / len(comp) * DX
                    cy = sum(i // NC for i in comp) / len(comp) * DY
                    svg.append(f'<text x="{cx:.1f}" y="{cy:.1f}" fill="#ffff00" font-size="13">E{eid}.{ci+1}</text>')
        # Hide labels in the dense full view; a separate review SVG carries them.
        svg.append("</svg>")
        (ROOT / f"slice001_{name.rsplit('.', 1)[0]}_reconstructed.svg").write_text("\n".join(svg))
        (ROOT / f"slice001_{name.rsplit('.', 1)[0]}_nolabel.svg").write_text("\n".join(svg))
    print(f"reconstructed_edges={len(confidence)} overlays=4")

if __name__ == "__main__":
    main()
