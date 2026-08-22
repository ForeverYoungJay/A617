#!/usr/bin/env bash
set -euo pipefail

ROOT=${A617_SERVER_ROOT:-$HOME/a617}
RAW="$ROOT/raw"
OUT="$ROOT/processed"
CODE="$ROOT/code"

python3 "$CODE/a617_stage1.py" --root "$RAW" --out "$OUT"
test -s "$OUT/slice_manifest.csv"
python3 - "$OUT/slice_manifest.csv" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
assert rows and all(int(r['points']) > 0 and int(r['eds_channels']) > 0 for r in rows)
print(f"validated slices={len(rows)}")
PY

# Raw instrument files are disposable only after the processed contract passes.
rm -rf -- "$RAW"
printf 'removed raw dataset; processed data remains in %s\n' "$OUT"
