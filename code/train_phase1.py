#!/usr/bin/env python3
"""Benchmark entry point; dependencies are intentionally explicit."""
from __future__ import annotations
import argparse, importlib.util
from pathlib import Path

def main():
    p = argparse.ArgumentParser(); p.add_argument("--data", type=Path, required=True); p.add_argument("--out", type=Path, required=True); a = p.parse_args()
    missing = [x for x in ("numpy", "sklearn", "torch") if importlib.util.find_spec(x) is None]
    if missing: raise SystemExit("install before training: " + ", ".join(missing))
    raise SystemExit("phase1 contract ready; implement Ridge/RF/MLP/GINE trainers against the CSV schema")

if __name__ == "__main__": main()
