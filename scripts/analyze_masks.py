#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from mask_analysis_common import (
    BACKBONE_DEFAULTS,
    bits_to_labels,
    compare_csv_to_summary,
    escape_summary,
    get_mask_counts,
    get_test_labels,
    labels_to_mask,
    load_hyb_run,
    minimal_hitting_sets,
    pair_witness_rows,
)
from verifier_mask_level import build_report


def parse_backbone_arg(value: str) -> tuple[str, str]:
    normalized = value.replace("{", "").replace("}", "").replace(" ", "")
    parts = [part for part in normalized.replace(";", ",").split(",") if part]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(f"Backbone must contain exactly two primes: {value!r}")
    return parts[0], parts[1]


def write_pair_witness_csv(path: Path, pair_rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "cover_prime_1",
                "cover_prime_2",
                "cover_mask",
                "passes",
                "witness_mask",
                "witness_count",
                "witness_kill_primes",
            ],
        )
        writer.writeheader()
        for row in pair_rows:
            writer.writerow({
                "cover_prime_1": row["cover_labels"][0],
                "cover_prime_2": row["cover_labels"][1],
                "cover_mask": row["cover_mask"],
                "passes": "true" if row["passes"] else "false",
                "witness_mask": "" if row["witness_mask"] is None else row["witness_mask"],
                "witness_count": row["witness_count"],
                "witness_kill_primes": ";".join(row["witness_kill_primes"]),
            })


def write_summary_markdown(
    path: Path,
    input_path: Path,
    run: dict,
    labels: list[str],
    mask_counts: dict[int, int],
    minimal_cover_size: int,
    minimal_covers: list[int],
    backbone_summaries: list[dict],
    csv_comparison: dict | None,
    checks: list[tuple[str, bool, str]],
) -> None:
    lines = [
        "# Mask-Level Anomaly Summary",
        "",
        f"- Input: `{input_path}`",
        f"- Origin: `{run.get('_origin', input_path)}`",
        f"- Completed: `{run.get('_completed', False)}`",
        f"- near_count: `{run.get('near_count')}`",
        f"- unique_mask_count: `{len(mask_counts)}`",
        f"- minimal_cover_size: `{minimal_cover_size}`",
        f"- minimal_covers_sample: `{'; '.join(','.join(bits_to_labels(mask, labels)) for mask in minimal_covers[:10])}`",
        "",
        "## Checks",
        "",
    ]
    for name, passed, detail in checks:
        lines.append(f"- {'PASS' if passed else 'FAIL'}: {name} ({detail})")

    lines.extend([
        "",
        "## Backbone Escapes",
        "",
    ])
    for summary in backbone_summaries:
        left, right = summary["backbone"]
        lines.append(f"### {{{left},{right}}}")
        lines.append("")
        lines.append(f"- Escaping near-candidates: `{summary['near_count']}`")
        lines.append(f"- Escaping unique masks: `{summary['unique_mask_count']}`")
        lines.append(f"- Overlap masks with empty kill set: `{summary['overlap_masks']}`")
        lines.append("- Third-prime mop-up distribution:")
        third_prime_counts = summary["third_prime_counts"]
        if third_prime_counts:
            for prime, count in third_prime_counts.items():
                lines.append(f"  - `{prime}`: `{count}`")
        else:
            lines.append("  - none")
        lines.append("- Escaping masks:")
        for row in summary["escaping_rows"]:
            lines.append(
                f"  - mask `{row['mask']}` count `{row['count']}` kills `{','.join(row['kill_primes']) or 'none'}`"
            )
        lines.append("")

    if len(backbone_summaries) >= 2:
        first_masks = {row["mask"] for row in backbone_summaries[0]["escaping_rows"]}
        second_masks = {row["mask"] for row in backbone_summaries[1]["escaping_rows"]}
        overlap = sorted(first_masks & second_masks)
        lines.append("## Backbone Overlap")
        lines.append("")
        lines.append(f"- Shared escaping masks between the first two backbones: `{overlap}`")
        lines.append("")

    if csv_comparison is not None:
        lines.append("## CSV Cross-Check")
        lines.append("")
        lines.append(f"- CSV rows: `{csv_comparison['row_count']}`")
        lines.append(f"- Expected rows from JSON: `{csv_comparison['expected_row_count']}`")
        if csv_comparison["mismatches"]:
            lines.append("- Mismatches:")
            for mismatch in csv_comparison["mismatches"]:
                lines.append(f"  - {mismatch}")
        else:
            lines.append("- CSV matches the JSON-derived escaping masks.")
        lines.append("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze mask-level anomaly behavior from a hyb JSON or ZIP artifact.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/input/erdos287-cover-hyb-cpp-N65-P65-S2-0.zip"),
        help="Path to the main N65 JSON or ZIP artifact",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("data/input/N65_backbone_escaping_masks.csv"),
        help="Optional helper CSV to cross-check",
    )
    parser.add_argument(
        "--backbone",
        action="append",
        type=parse_backbone_arg,
        default=None,
        help="Repeatable backbone pair like 2,31",
    )
    parser.add_argument("--expected-near-count", type=int, default=None)
    parser.add_argument("--expected-unique-mask-count", type=int, default=None)
    parser.add_argument("--expected-minimal-cover-size", type=int, default=None)
    parser.add_argument("--expected-escape-count", action="append", default=None, metavar="PAIR=COUNT")
    parser.add_argument(
        "--expect-all-two-prime-covers-fail",
        action="store_true",
        help="Require every 2-prime cover to fail with a witness mask",
    )
    parser.add_argument("--out-summary", type=Path, default=Path("results/mask_anomaly_summary.md"))
    parser.add_argument("--out-witnesses", type=Path, default=Path("results/minimality_witnesses.csv"))
    parser.add_argument("--out-log", type=Path, default=Path("results/mask_verifier_log.txt"))
    args = parser.parse_args()

    backbones = args.backbone or list(BACKBONE_DEFAULTS)
    expected_escape_counts: dict[tuple[str, str], int] = {}
    for item in args.expected_escape_count or []:
        pair_text, count_text = item.split("=", 1)
        expected_escape_counts[parse_backbone_arg(pair_text)] = int(count_text)

    run = load_hyb_run(args.input)
    labels = get_test_labels(run)
    mask_counts = get_mask_counts(run)
    minimal_cover_size, minimal_covers = minimal_hitting_sets(mask_counts, len(labels))
    pair_rows = pair_witness_rows(mask_counts, labels)
    backbone_summaries = [escape_summary(mask_counts, labels, backbone) for backbone in backbones]

    checks: list[tuple[str, bool, str]] = []
    checks.append(("artifact completed", bool(run.get("_completed", False)), f"completed={run.get('_completed', False)}"))
    if args.expected_near_count is not None:
        checks.append((
            "near_count",
            int(run.get("near_count", -1)) == args.expected_near_count,
            f"expected={args.expected_near_count} actual={run.get('near_count')}",
        ))
    if args.expected_unique_mask_count is not None:
        checks.append((
            "unique_mask_count",
            len(mask_counts) == args.expected_unique_mask_count,
            f"expected={args.expected_unique_mask_count} actual={len(mask_counts)}",
        ))
    if args.expected_minimal_cover_size is not None:
        checks.append((
            "minimal_cover_size",
            minimal_cover_size == args.expected_minimal_cover_size,
            f"expected={args.expected_minimal_cover_size} actual={minimal_cover_size}",
        ))
    for summary in backbone_summaries:
        backbone = tuple(summary["backbone"])
        if backbone in expected_escape_counts:
            checks.append((
                f"escape_count {{{backbone[0]},{backbone[1]}}}",
                summary["near_count"] == expected_escape_counts[backbone],
                f"expected={expected_escape_counts[backbone]} actual={summary['near_count']}",
            ))

    all_pairs_fail = all(not row["passes"] for row in pair_rows)
    if args.expect_all_two_prime_covers_fail:
        checks.append(("all 2-prime covers fail", all_pairs_fail, f"pair_count={len(pair_rows)}"))

    csv_comparison = None
    if args.csv.exists():
        csv_comparison = compare_csv_to_summary(args.csv, backbone_summaries, labels)
        checks.append((
            "helper CSV matches JSON",
            not csv_comparison["mismatches"],
            f"mismatches={len(csv_comparison['mismatches'])}",
        ))
    else:
        checks.append(("helper CSV present", True, f"skipped missing={args.csv}"))

    write_pair_witness_csv(args.out_witnesses, pair_rows)
    args.out_log.parent.mkdir(parents=True, exist_ok=True)
    args.out_log.write_text(build_report(args.input), encoding="utf-8")
    write_summary_markdown(
        args.out_summary,
        args.input,
        run,
        labels,
        mask_counts,
        minimal_cover_size,
        minimal_covers,
        backbone_summaries,
        csv_comparison,
        checks,
    )

    for name, passed, detail in checks:
        print(f"{'PASS' if passed else 'FAIL'}: {name} ({detail})")

    if not all(passed for _, passed, _ in checks):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
