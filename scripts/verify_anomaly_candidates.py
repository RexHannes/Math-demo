#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter
from fractions import Fraction
from pathlib import Path


def primes_upto(limit: int) -> list[int]:
    if limit < 2:
        return []
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for value in range(2, int(limit**0.5) + 1):
        if not sieve[value]:
            continue
        for multiple in range(value * value, limit + 1, value):
            sieve[multiple] = False
    return [value for value in range(2, limit + 1) if sieve[value]]


def parse_expected_count(value: str) -> tuple[str, int]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("Expected count must look like 2+31=NN")
    key, count_text = value.split("=", 1)
    return key.strip(), int(count_text)


def parse_denominators(text: str) -> list[int]:
    separator = ";" if ";" in text else ","
    return [int(value) for value in text.split(separator) if value.strip()]


def row_backbones(row: dict[str, str]) -> list[str]:
    return [value for value in row.get("escaped_backbone", "").split("|") if value]


def main() -> None:
    parser = argparse.ArgumentParser(description="Independently verify merged anomaly candidate rows.")
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--N", type=int, default=None)
    parser.add_argument("--P", type=int, default=65)
    parser.add_argument("--low", default="999/1000")
    parser.add_argument("--high", default="1001/1000")
    parser.add_argument("--out-report", type=Path, default=Path("results/anomaly_candidate_verification.md"))
    parser.add_argument("--expected-total", type=int, default=None)
    parser.add_argument("--expected-backbone-count", action="append", type=parse_expected_count, default=None)
    args = parser.parse_args()

    labels = [str(prime) for prime in primes_upto(args.P)]
    prime_to_bit = {label: index for index, label in enumerate(labels)}
    low = Fraction(args.low)
    high = Fraction(args.high)
    expected_counts = dict(args.expected_backbone_count or [])

    rows: list[dict[str, str]] = []
    errors: list[str] = []
    backbone_counts: Counter[str] = Counter()

    with args.csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row_number, row in enumerate(reader, start=2):
            rows.append(row)
            denominators = parse_denominators(row.get("denominators", ""))
            if not denominators:
                errors.append(f"row {row_number}: no denominators")
                continue

            exact = sum((Fraction(1, denominator) for denominator in denominators), Fraction(0, 1))
            if not low <= exact <= high:
                errors.append(f"row {row_number}: exact sum {exact} outside [{low}, {high}]")

            try:
                kill_mask = int(row.get("kill_mask", ""))
            except ValueError:
                errors.append(f"row {row_number}: invalid kill_mask `{row.get('kill_mask', '')}`")
                continue

            for backbone in row_backbones(row):
                backbone_counts[backbone] += 1
                for prime_label in backbone.split("+"):
                    bit = prime_to_bit.get(prime_label)
                    if bit is None:
                        errors.append(f"row {row_number}: backbone prime {prime_label} is not in P={args.P}")
                    elif kill_mask & (1 << bit):
                        errors.append(f"row {row_number}: row is killed by escaped backbone prime {prime_label}")

    if args.expected_total is not None and len(rows) != args.expected_total:
        errors.append(f"expected {args.expected_total} rows, found {len(rows)}")
    for backbone, expected_count in sorted(expected_counts.items()):
        actual = backbone_counts[backbone]
        if actual != expected_count:
            errors.append(f"expected {expected_count} rows for {backbone}, found {actual}")

    args.out_report.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"# N{args.N if args.N is not None else 'unknown'} Anomaly Candidate Verification",
        "",
        f"- N: `{args.N if args.N is not None else 'unknown'}`",
        f"- P: `{args.P}`",
        f"- window: `[{args.low}, {args.high}]`",
        f"- csv: `{args.csv_path}`",
        f"- rows: `{len(rows)}`",
        f"- expected_rows: `{args.expected_total if args.expected_total is not None else 'not enforced'}`",
        f"- errors: `{len(errors)}`",
        "",
        "## Backbone Counts",
        "",
    ]
    if expected_counts:
        for backbone, expected_count in sorted(expected_counts.items()):
            lines.append(f"- `{backbone}`: `{backbone_counts[backbone]}` / `{expected_count}`")
    elif backbone_counts:
        for backbone, count in sorted(backbone_counts.items()):
            lines.append(f"- `{backbone}`: `{count}`")
    else:
        lines.append("- none")
    lines.extend(["", "## Checks", ""])
    lines.append(f"- exact sums in window: `{'yes' if not any('exact sum' in error for error in errors) else 'no'}`")
    lines.append(f"- escaped backbone bits absent from kill masks: `{'yes' if not any('killed by escaped' in error for error in errors) else 'no'}`")
    if errors:
        lines.extend(["", "## Errors", ""])
        lines.extend(f"- {error}" for error in errors[:100])
        if len(errors) > 100:
            lines.append(f"- ... {len(errors) - 100} more")
    args.out_report.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Verified {len(rows)} rows from {args.csv_path}")
    print(f"Errors: {len(errors)}")
    if expected_counts:
        for backbone, expected_count in sorted(expected_counts.items()):
            print(f"{backbone}: {backbone_counts[backbone]}/{expected_count}")
    else:
        for backbone, count in sorted(backbone_counts.items()):
            print(f"{backbone}: {count}")
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
