#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from mask_analysis_common import (
    bits_to_labels,
    get_mask_counts,
    get_test_labels,
    load_hyb_run,
    minimal_hitting_sets,
    pair_witness_rows,
)


def build_report(input_path: Path) -> str:
    run = load_hyb_run(input_path)
    labels = get_test_labels(run)
    mask_counts = get_mask_counts(run)
    min_size, solutions = minimal_hitting_sets(mask_counts, len(labels))
    pair_rows = pair_witness_rows(mask_counts, labels)
    failing_pairs = [row for row in pair_rows if not row["passes"]]
    passing_pairs = [row for row in pair_rows if row["passes"]]

    lines = [
        f"input={input_path}",
        f"origin={run.get('_origin', input_path)}",
        f"completed={run.get('_completed', False)}",
        f"near_count={run.get('near_count')}",
        f"unique_mask_count={len(mask_counts)}",
        f"minimal_cover_size={min_size}",
        "minimal_covers_sample=" + "; ".join(",".join(bits_to_labels(mask, labels)) for mask in solutions[:10]),
        f"pair_count={len(pair_rows)}",
        f"pair_failures={len(failing_pairs)}",
        f"pair_passes={len(passing_pairs)}",
    ]
    if passing_pairs:
        lines.append(
            "unexpected_pair_passes="
            + "; ".join(",".join(row["cover_labels"]) for row in passing_pairs[:10])
        )
    else:
        lines.append("all_two_prime_covers_fail=true")
    if failing_pairs:
        sample = failing_pairs[:10]
        lines.append(
            "sample_witnesses="
            + "; ".join(
                f"{row['cover_labels'][0]},{row['cover_labels'][1]} -> {row['witness_mask']} ({','.join(row['witness_kill_primes'])})"
                for row in sample
            )
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Independently verify mask-level minimal cover size from unique_masks.")
    parser.add_argument("input", type=Path, help="Path to a hyb result JSON or ZIP artifact")
    parser.add_argument("--out", type=Path, default=None, help="Optional text log path")
    args = parser.parse_args()

    report = build_report(args.input)
    print(report, end="")
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
