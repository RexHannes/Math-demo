#!/usr/bin/env python3
import csv
import sys
from pathlib import Path

if len(sys.argv) < 2:
    print("Usage: python scripts/summarize_modulus_probe.py <artifact_or_results_dir>")
    raise SystemExit(2)

root = Path(sys.argv[1])
files = sorted(root.rglob("*_summary_by_M.csv"))
if not files:
    print(f"No *_summary_by_M.csv files found under {root}")
    raise SystemExit(1)

rows_by_m = {}
def parse_int(value):
    return int(value) if value not in ("", None) else -1

for f in files:
    try:
        with f.open(newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                row["M"] = parse_int(row["M"])
                row["mstar"] = parse_int(row["mstar"])
                row["total"] = parse_int(row["total"])
                row["blocked"] = parse_int(row["blocked"])
                row["reachable"] = parse_int(row["reachable"])
                row["first_reachable_m"] = parse_int(row["first_reachable_m"])
                row["last_reachable_m"] = parse_int(row["last_reachable_m"])
                row["full_block"] = row["blocked"] == row["total"]
                row["partial_block"] = row["blocked"] > 0 and not row["full_block"]
                rows_by_m[row["M"]] = row
    except Exception as exc:
        print(f"Skipping {f}: {exc}")

rows = [rows_by_m[m] for m in sorted(rows_by_m)]
print("Overall:")
print({
    "M_min": rows[0]["M"],
    "M_max": rows[-1]["M"],
    "rows": len(rows),
    "full_block_M_count": sum(1 for row in rows if row["full_block"]),
    "partial_block_M_count": sum(1 for row in rows if row["partial_block"]),
})

runs = []
start = prev = None
for row in rows:
    m_value = row["M"]
    if row["full_block"]:
        if start is None:
            start = prev = m_value
        elif m_value == prev + 1:
            prev = m_value
        else:
            runs.append((start, prev))
            start = prev = m_value
    else:
        if start is not None:
            runs.append((start, prev))
            start = prev = None
if start is not None:
    runs.append((start, prev))

print("\nFull-block runs:")
for start_value, end_value in runs:
    print(f"  {start_value}-{end_value} ({end_value-start_value+1} values)")

print("\nRows around regime boundaries:")
interesting = [
    row for row in rows
    if 230 <= row["M"] <= 250 or 355 <= row["M"] <= 370
]
header = ["M", "mstar", "total", "blocked", "reachable", "first_reachable_m", "last_reachable_m"]
print(" ".join(f"{name:>17}" for name in header))
for row in interesting:
    print(" ".join(f"{str(row[name]):>17}" for name in header))

out = root / "merged_modulus_probe_summary.csv"
with out.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=[
        "M", "mstar", "total", "blocked", "reachable",
        "first_reachable_m", "last_reachable_m",
        "example_blocked_nz", "example_reachable_nz",
        "full_block", "partial_block"
    ])
    writer.writeheader()
    writer.writerows(rows)
print(f"\nWrote {out}")
