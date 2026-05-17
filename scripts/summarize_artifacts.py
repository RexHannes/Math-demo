#!/usr/bin/env python3
"""Merge GitHub Actions artifact CSV/JSON outputs for the Erdős interval diagnostic.
Usage:
  python scripts/summarize_artifacts.py path/to/unzipped_artifacts
or from repo root after downloading/extracting artifacts into ./artifacts:
  python scripts/summarize_artifacts.py artifacts
"""
from __future__ import annotations
import csv, json, sys
from pathlib import Path
from collections import Counter, defaultdict

root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("artifacts")
if not root.exists():
    raise SystemExit(f"Not found: {root}")

summary_files = sorted(root.rglob("*_summary_by_M.csv"))
record_files = sorted(root.rglob("*_records.csv"))
report_files = sorted(root.rglob("*_report.json"))
print(f"summary files: {len(summary_files)}")
print(f"record files:  {len(record_files)}")
print(f"report files:  {len(report_files)}")

by_M = {}
overall = Counter()
for f in summary_files:
    with f.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            M = int(row["M"])
            # One file per M should have one row, but combine safely.
            acc = by_M.setdefault(M, Counter())
            for k, v in row.items():
                if k == "M":
                    continue
                acc[k] += int(v or 0)

for M, c in by_M.items():
    for k, v in c.items():
        if k != "total":
            overall[k] += v
    overall["total"] += c.get("total", 0)

print("\nOverall counts")
for k in ["total", "direct_prime", "mixed_or_large_modulus", "dense_forced_holes", "linked_prime_pair", "forced_modular", "endpoint_forced", "no_state_forced", "survivor"]:
    print(f"{k:26s} {overall[k]}")

print("\nBy M")
for M in sorted(by_M):
    c = by_M[M]
    print(f"M={M:3d} total={c['total']:4d} prime={c['direct_prime']:4d} mixed={c['mixed_or_large_modulus']:4d} dense={c['dense_forced_holes']:4d} linked={c['linked_prime_pair']:4d} fmod={c['forced_modular']:4d} survivor={c['survivor']:4d}")

# Top linked pairs and signatures, if any reports contain them.
pairs = Counter()
sigs = Counter()
for f in report_files:
    try:
        data = json.loads(f.read_text())
    except Exception:
        continue
    pairs.update(data.get("top_linked_pairs", {}))
    sigs.update(data.get("top_signatures", {}))

print("\nTop linked pairs")
for k, v in pairs.most_common(20):
    print(f"{k:12s} {v}")

print("\nTop linked signatures")
for k, v in sigs.most_common(20):
    print(f"{k:45s} {v}")

out = root / "merged_summary_by_M.csv"
fieldnames = ["M", "total", "direct_prime", "mixed_or_large_modulus", "dense_forced_holes", "linked_prime_pair", "forced_modular", "endpoint_forced", "no_state_forced", "survivor"]
with out.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames)
    writer.writeheader()
    for M in sorted(by_M):
        row = {"M": M}
        row.update({k: by_M[M].get(k, 0) for k in fieldnames if k != "M"})
        writer.writerow(row)
print(f"\nWrote {out}")
