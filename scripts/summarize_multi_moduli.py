#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path


def find_windows(rows, min_len):
    windows = []
    start = prev = None
    for row in rows:
        m_value = row["M"]
        full = row["full_block"]
        if full and start is None:
            start = prev = m_value
        elif full and m_value == prev + 1:
            prev = m_value
        elif full:
            if prev - start + 1 >= min_len:
                windows.append((start, prev, prev - start + 1))
            start = prev = m_value
        elif start is not None:
            if prev - start + 1 >= min_len:
                windows.append((start, prev, prev - start + 1))
            start = prev = None
    if start is not None and prev - start + 1 >= min_len:
        windows.append((start, prev, prev - start + 1))
    return windows


def parse_summary_file(path):
    match = re.search(r"(?:mod|probe)(\d+)(?:_M|_)", path.name)
    if not match:
        return []
    modulus = int(match.group(1))
    rows = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            parsed = {
                "modulus": modulus,
                "M": int(row["M"]),
                "mstar": int(row["mstar"]),
                "total": int(row["total"]),
                "blocked": int(row["blocked"]),
                "reachable": int(row["reachable"]),
                "first_reachable_m": int(row["first_reachable_m"]) if row["first_reachable_m"] else -1,
                "last_reachable_m": int(row["last_reachable_m"]) if row["last_reachable_m"] else -1,
                "example_blocked_nz": row.get("example_blocked_nz", ""),
                "example_reachable_nz": row.get("example_reachable_nz", "")
            }
            parsed["full_block"] = parsed["blocked"] == parsed["total"]
            parsed["partial_block"] = parsed["blocked"] > 0 and parsed["blocked"] < parsed["total"]
            parsed["block_rate"] = (parsed["blocked"] / parsed["total"]) if parsed["total"] else 0.0
            rows.append(parsed)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("artifacts_dir", help="Directory containing unzipped GitHub artifact CSVs")
    parser.add_argument("--min-len", type=int, default=20, help="Minimum full-block window length to report")
    parser.add_argument("--out-prefix", default="multi_moduli_summary")
    args = parser.parse_args()

    root = Path(args.artifacts_dir)
    csvs = sorted(root.rglob("*_summary_by_M.csv"))
    if not csvs:
        raise SystemExit(f"No *_summary_by_M.csv found under {root}")

    by_key = {}
    for path in csvs:
        for row in parse_summary_file(path):
            by_key[(row["modulus"], row["M"])] = row

    if not by_key:
        raise SystemExit("No parseable modulus summary files found")

    all_rows = [by_key[key] for key in sorted(by_key)]

    all_by_m_path = Path(f"{args.out_prefix}_all_by_M.csv")
    with all_by_m_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "modulus", "M", "mstar", "total", "blocked", "reachable",
            "first_reachable_m", "last_reachable_m",
            "example_blocked_nz", "example_reachable_nz",
            "full_block", "partial_block", "block_rate"
        ])
        writer.writeheader()
        writer.writerows(all_rows)

    summary_rows = []
    window_rows = []
    grouped = {}
    for row in all_rows:
        grouped.setdefault(row["modulus"], []).append(row)

    for modulus, rows in sorted(grouped.items()):
        rows.sort(key=lambda item: item["M"])
        summary_rows.append({
            "modulus": modulus,
            "M_min": rows[0]["M"],
            "M_max": rows[-1]["M"],
            "num_M": len(rows),
            "total_intervals": sum(row["total"] for row in rows),
            "blocked_intervals": sum(row["blocked"] for row in rows),
            "reachable_intervals": sum(row["reachable"] for row in rows),
            "full_block_M_count": sum(1 for row in rows if row["full_block"]),
            "partial_block_M_count": sum(1 for row in rows if row["partial_block"]),
            "max_block_rate": max(row["block_rate"] for row in rows)
        })
        for start, end, length in find_windows(rows, args.min_len):
            window_rows.append({
                "modulus": modulus,
                "start_M": start,
                "end_M": end,
                "length": length
            })

    summary_rows.sort(key=lambda row: (-row["full_block_M_count"], -row["blocked_intervals"], row["modulus"]))
    window_rows.sort(key=lambda row: (-row["length"], row["modulus"], row["start_M"]))

    by_modulus_path = Path(f"{args.out_prefix}_by_modulus.csv")
    with by_modulus_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "modulus", "M_min", "M_max", "num_M", "total_intervals",
            "blocked_intervals", "reachable_intervals", "full_block_M_count",
            "partial_block_M_count", "max_block_rate"
        ])
        writer.writeheader()
        writer.writerows(summary_rows)

    windows_path = Path(f"{args.out_prefix}_windows_minlen{args.min_len}.csv")
    with windows_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["modulus", "start_M", "end_M", "length"])
        writer.writeheader()
        writer.writerows(window_rows)

    print("\n=== By modulus ===")
    for row in summary_rows:
        print(row)
    print(f"\n=== Full-block windows length >= {args.min_len} ===")
    for row in window_rows:
        print(row)


if __name__ == "__main__":
    main()
