#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


def parse_expected_count(value: str) -> tuple[str, int]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("Expected count must look like 2+31=NN")
    key, count_text = value.split("=", 1)
    return key.strip(), int(count_text)


def parse_chunk_label(path: Path) -> str:
    stem = path.stem
    marker = "_S"
    if marker in stem:
        return "S" + stem.split(marker, 1)[1]
    return stem


def load_rows(paths: list[Path]) -> tuple[list[dict[str, str]], list[str]]:
    rows: list[dict[str, str]] = []
    fieldnames: list[str] = []
    for path in sorted(paths):
        chunk_label = parse_chunk_label(path)
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames and not fieldnames:
                fieldnames = list(reader.fieldnames) + ["source_chunk", "source_csv"]
            for index, row in enumerate(reader, start=1):
                merged = dict(row)
                merged["source_chunk"] = chunk_label
                merged["source_csv"] = str(path)
                original_id = (row.get("candidate_id") or "").strip() or str(index)
                merged["candidate_id"] = f"{chunk_label}:{original_id}"
                rows.append(merged)
    return rows, fieldnames


def load_statuses(paths: list[Path]) -> list[dict[str, object]]:
    statuses: list[dict[str, object]] = []
    for path in sorted(paths):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            payload = {
                "chunk_id": parse_chunk_label(path),
                "status": "invalid",
                "error": str(exc),
            }
        payload["source_status"] = str(path)
        statuses.append(payload)
    return statuses


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        if not fieldnames and rows:
            fieldnames = list(rows[0].keys())
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        if rows:
            writer.writerows(rows)


def parse_int_list(text: str) -> list[int]:
    if not text:
        return []
    return [int(value) for value in text.split(";") if value.strip()]


def row_backbones(row: dict[str, str]) -> list[str]:
    return [value for value in row.get("escaped_backbone", "").split("|") if value]


def row_denominators(row: dict[str, str]) -> list[int]:
    return parse_int_list(row.get("denominators", ""))


def classify_n65_template(row: dict[str, str], denominators: list[int]) -> str:
    denominator_set = set(denominators)
    if denominators and max(denominators) > 65:
        return ""
    backbones = set(row_backbones(row))
    if "2+31" in backbones and 16 not in denominator_set:
        if {15, 17, 19}.issubset(denominator_set):
            return "N65-T1 2+31 miss16 core15_17_19"
        return "N65-T2 2+31 miss16 relaxed"
    if "19+37" in backbones and 19 not in denominator_set:
        if {18, 20, 22}.issubset(denominator_set):
            return "N65-T3 19+37 miss19 core18_20_22"
        return "N65-T4 19+37 miss19 relaxed"
    return ""


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def counter_lines(counter: Counter, label: str, limit: int | None = None) -> list[str]:
    lines = [f"## {label}", ""]
    items = counter.most_common(limit)
    if not items:
        lines.append("- none")
    for key, count in items:
        lines.append(f"- `{key}`: `{count}`")
    lines.append("")
    return lines


def write_summary(
    path: Path,
    rows: list[dict[str, str]],
    statuses: list[dict[str, object]],
    expected_counts: dict[str, int],
    expected_total: int | None,
    expected_chunk_count: int | None,
    denominator_max: int,
    run_N: int | None,
    run_P: int | None,
    low: str | None,
    high: str | None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    backbone_counts: Counter[str] = Counter()
    unique_masks_by_backbone: defaultdict[str, set[str]] = defaultdict(set)
    chunk_counts: Counter[str] = Counter()
    length_counts: Counter[int] = Counter()
    interval_counts: Counter[str] = Counter()
    gap_counts: Counter[str] = Counter()
    denominator_counts: defaultdict[str, Counter[int]] = defaultdict(Counter)
    missing_denominator_counts: defaultdict[str, Counter[int]] = defaultdict(Counter)
    structural_flags = Counter()
    third_prime_counts: defaultdict[str, Counter[str]] = defaultdict(Counter)
    rows_by_backbone: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    max_denominator_counts: Counter[int] = Counter()
    rows_with_max_above_65: list[dict[str, str]] = []
    template_counts: Counter[str] = Counter()
    unclassified_rows: list[dict[str, str]] = []
    universe = set(range(2, denominator_max + 1))

    for row in rows:
        chunk_counts[row["source_chunk"]] += 1
        try:
            length_counts[int(row["length"])] += 1
        except Exception:
            pass

        escaped_backbones = row_backbones(row)
        denominators = row_denominators(row)
        denominator_set = set(denominators)
        if denominators:
            max_denominator = max(denominators)
            max_denominator_counts[max_denominator] += 1
            if max_denominator > 65:
                rows_with_max_above_65.append(row)
        template = classify_n65_template(row, denominators)
        if template:
            template_counts[template] += 1
        else:
            unclassified_rows.append(row)
        interval = f"{row.get('min_denominator', '')}-{row.get('max_denominator', '')}"
        if interval != "-":
            interval_counts[interval] += 1
        gap_pattern = row.get("gap_pattern", "")
        if gap_pattern:
            gap_counts[gap_pattern] += 1

        for backbone in escaped_backbones:
            backbone_counts[backbone] += 1
            rows_by_backbone[backbone].append(row)
            if row.get("kill_mask"):
                unique_masks_by_backbone[backbone].add(row["kill_mask"])
            denominator_counts[backbone].update(denominators)
            missing_denominator_counts[backbone].update(sorted(universe - denominator_set))

        if row.get("contains_37_38") == "true":
            structural_flags["contains_37_38"] += 1
        if row.get("contains_61_62") == "true":
            structural_flags["contains_61_62"] += 1
        if row.get("contains_adjacent_factor_pair") == "true":
            structural_flags["contains_adjacent_factor_pair"] += 1

        kill_primes = [value for value in row.get("kill_primes", "").split(";") if value]
        for backbone in escaped_backbones:
            backbone_primes = set(backbone.split("+"))
            for prime in kill_primes:
                if prime not in backbone_primes:
                    third_prime_counts[backbone][prime] += 1

    row_count_complete = expected_total is None or len(rows) == expected_total
    backbone_counts_complete = all(backbone_counts[key] == value for key, value in expected_counts.items())
    complete_status_count = sum(1 for status in statuses if status.get("status") == "complete")
    incomplete_statuses = [status for status in statuses if status.get("status") != "complete"]
    chunk_count_complete = expected_chunk_count is None or len(statuses) == expected_chunk_count
    chunk_status_complete = bool(statuses) and chunk_count_complete and not incomplete_statuses
    commit_shas = {str(status.get("commit_sha", "")) for status in statuses if status.get("commit_sha")}
    source_hashes = {str(status.get("cpp_sha256", "")) for status in statuses if status.get("cpp_sha256")}
    compile_flags = {str(status.get("compile_flags", "")) for status in statuses if status.get("compile_flags")}
    build_metadata_consistent = len(commit_shas) <= 1 and len(source_hashes) <= 1 and len(compile_flags) <= 1
    complete = chunk_status_complete and row_count_complete and backbone_counts_complete and build_metadata_consistent

    lines = [
        f"# N{run_N if run_N is not None else 'unknown'} Anomaly Candidate Summary",
        "",
        f"- N: `{run_N if run_N is not None else 'unknown'}`",
        f"- P: `{run_P if run_P is not None else 'unknown'}`",
        f"- window: `[{low if low is not None else 'unknown'}, {high if high is not None else 'unknown'}]`",
        f"- total anomaly rows found: `{len(rows)}`",
        f"- expected_rows: `{expected_total if expected_total is not None else 'not enforced'}`",
        f"- complete: `{yes_no(complete)}`",
        f"- complete_chunk_status: `{yes_no(chunk_status_complete)}`",
        f"- chunk_status_files: `{len(statuses)}`",
        f"- expected_chunk_status_files: `{expected_chunk_count if expected_chunk_count is not None else 'not enforced'}`",
        f"- complete_chunk_status_files: `{complete_status_count}`",
        f"- distinct_commit_shas: `{len(commit_shas)}`",
        f"- distinct_cpp_sha256: `{len(source_hashes)}`",
        f"- distinct_compile_flags: `{len(compile_flags)}`",
        "",
        "## Completeness",
        "",
        f"- total rows found: `{len(rows)}`",
    ]
    if expected_total is not None:
        lines.append(f"- total rows expected: `{expected_total}`")
    for backbone, expected_count in sorted(expected_counts.items()):
        lines.append(f"- `{backbone}` rows: `{backbone_counts[backbone]}` / `{expected_count}`")

    lines.extend([
        "",
        "## Chunk Status Sidecars",
        "",
    ])
    if statuses:
        for status in sorted(statuses, key=lambda item: str(item.get("chunk_id", ""))):
            chunk_id = status.get("chunk_id", "unknown")
            state = status.get("status", "unknown")
            candidates = status.get("candidates_checked", "unknown")
            exceptions = status.get("exceptions_found", "unknown")
            lines.append(
                f"- `{chunk_id}`: status=`{state}` candidates_checked=`{candidates}` exceptions_found=`{exceptions}`"
            )
    else:
        lines.append("- none")

    if incomplete_statuses:
        lines.extend([
            "",
            "## Incomplete Or Invalid Chunks",
            "",
        ])
        for status in sorted(incomplete_statuses, key=lambda item: str(item.get("chunk_id", ""))):
            chunk_id = status.get("chunk_id", "unknown")
            state = status.get("status", "unknown")
            source = status.get("source_status", "")
            lines.append(f"- `{chunk_id}`: `{state}` source=`{source}`")

    lines.extend([
        "",
        "## All Source Chunks",
        "",
    ])
    source_chunks = sorted({str(status.get("chunk_id", "unknown")) for status in statuses} | set(chunk_counts))
    if source_chunks:
        for chunk in source_chunks:
            lines.append(f"- `{chunk}`")
    else:
        lines.append("- none")

    lines.extend([
        "",
        "## Row Count By Backbone",
        "",
    ])
    if backbone_counts:
        for backbone, count in sorted(backbone_counts.items()):
            lines.append(f"- `{backbone}`: `{count}`")
    else:
        lines.append("- none")

    lines.extend([
        "",
        "## Unique Kill Masks By Backbone",
        "",
    ])
    if unique_masks_by_backbone:
        for backbone, masks in sorted(unique_masks_by_backbone.items()):
            lines.append(f"- `{backbone}`: `{len(masks)}`")
    else:
        lines.append("- none")

    lines.extend([
        "",
        "## Chunk Counts",
        "",
    ])
    if chunk_counts:
        for chunk, count in sorted(chunk_counts.items()):
            lines.append(f"- `{chunk}`: `{count}`")
    else:
        lines.append("- none")

    lines.extend([
        "",
        "## Structural Flags",
        "",
        f"- contains_37_38: `{structural_flags['contains_37_38']}`",
        f"- contains_61_62: `{structural_flags['contains_61_62']}`",
        f"- contains_adjacent_factor_pair: `{structural_flags['contains_adjacent_factor_pair']}`",
        "",
        "## Length Distribution",
        "",
    ])
    if length_counts:
        for length, count in sorted(length_counts.items()):
            lines.append(f"- length `{length}`: `{count}`")
    else:
        lines.append("- none")

    lines.extend([
        "",
        "## Max Denominator Distribution",
        "",
    ])
    if max_denominator_counts:
        for denominator, count in sorted(max_denominator_counts.items()):
            lines.append(f"- `{denominator}`: `{count}`")
    else:
        lines.append("- none")

    lines.extend([
        "",
        "## Rows With Max Denominator Above 65",
        "",
    ])
    if rows_with_max_above_65:
        for row in rows_with_max_above_65[:50]:
            lines.append(
                f"- `{row['candidate_id']}` backbone=`{row.get('escaped_backbone', '')}` "
                f"max=`{row.get('max_denominator', '')}` denominators=`{row.get('denominators', '')}`"
            )
        if len(rows_with_max_above_65) > 50:
            lines.append(f"- ... `{len(rows_with_max_above_65) - 50}` more")
    else:
        lines.append("- none")

    lines.extend([
        "",
        "## N65 Template Classification Counts",
        "",
        "- `N65-T1 2+31 miss16 core15_17_19`: `{}`".format(template_counts["N65-T1 2+31 miss16 core15_17_19"]),
        "- `N65-T2 2+31 miss16 relaxed`: `{}`".format(template_counts["N65-T2 2+31 miss16 relaxed"]),
        "- `N65-T3 19+37 miss19 core18_20_22`: `{}`".format(template_counts["N65-T3 19+37 miss19 core18_20_22"]),
        "- `N65-T4 19+37 miss19 relaxed`: `{}`".format(template_counts["N65-T4 19+37 miss19 relaxed"]),
        "",
        "## N65 Template Definitions",
        "",
        "- `N65-T1`: escapes `2+31`, max denominator <= `65`, misses `16`, includes `15;17;19`",
        "- `N65-T2`: escapes `2+31`, max denominator <= `65`, misses `16`, but not the full `15;17;19` core",
        "- `N65-T3`: escapes `19+37`, max denominator <= `65`, misses `19`, includes `18;20;22`",
        "- `N65-T4`: escapes `19+37`, max denominator <= `65`, misses `19`, but not the full `18;20;22` core",
        "",
        "## Unclassified Rows",
        "",
    ])
    if unclassified_rows:
        for row in unclassified_rows[:50]:
            lines.append(
                f"- `{row['candidate_id']}` backbone=`{row.get('escaped_backbone', '')}` "
                f"max=`{row.get('max_denominator', '')}` denominators=`{row.get('denominators', '')}`"
            )
        if len(unclassified_rows) > 50:
            lines.append(f"- ... `{len(unclassified_rows) - 50}` more")
    else:
        lines.append("- none")
    lines.append("")

    lines.extend(counter_lines(interval_counts, "Interval Distribution"))
    lines.extend(counter_lines(gap_counts, "Repeated Gap Patterns", limit=30))

    lines.extend([
        "## Denominator Frequency",
        "",
    ])
    for backbone in sorted(rows_by_backbone):
        lines.append(f"### {backbone}")
        lines.append("")
        for denominator, count in sorted(denominator_counts[backbone].items()):
            lines.append(f"- `{denominator}`: `{count}`")
        lines.append("")

    lines.extend([
        "## Missing-Denominator Frequency",
        "",
    ])
    for backbone in sorted(rows_by_backbone):
        lines.append(f"### {backbone}")
        lines.append("")
        for denominator, count in sorted(missing_denominator_counts[backbone].items()):
            lines.append(f"- `{denominator}`: `{count}`")
        lines.append("")

    lines.extend([
        "## Universal Included And Missing Denominators",
        "",
    ])
    for backbone in sorted(rows_by_backbone):
        backbone_rows = rows_by_backbone[backbone]
        included_sets = [set(row_denominators(row)) for row in backbone_rows]
        universal_included = sorted(set.intersection(*included_sets)) if included_sets else []
        universal_missing = sorted(universe - set.union(*included_sets)) if included_sets else sorted(universe)
        lines.append(f"### {backbone}")
        lines.append("")
        lines.append(f"- universal included: `{';'.join(str(value) for value in universal_included) or 'none'}`")
        lines.append(f"- universal missing: `{';'.join(str(value) for value in universal_missing) or 'none'}`")
        if backbone == "2+31":
            misses_16 = all(16 not in row_denominators(row) for row in backbone_rows)
            lines.append(f"- every row misses 16: `{yes_no(misses_16)}`")
        if backbone == "19+37":
            misses_19 = all(19 not in row_denominators(row) for row in backbone_rows)
            lines.append(f"- every row misses 19: `{yes_no(misses_19)}`")
        lines.append("")

    lines.extend([
        "",
        "## Third-Prime Mop-Up Counts",
        "",
    ])
    if third_prime_counts:
        for backbone, counts in sorted(third_prime_counts.items()):
            lines.append(f"### {backbone}")
            lines.append("")
            for prime, count in sorted(counts.items(), key=lambda item: int(item[0])):
                lines.append(f"- `{prime}`: `{count}`")
            lines.append("")
    else:
        lines.append("- none")
        lines.append("")

    lines.extend([
        "## Sample Rows",
        "",
    ])
    for row in rows[:10]:
        lines.append(
            f"- `{row['candidate_id']}` backbone=`{row.get('escaped_backbone', '')}` "
            f"denominators=`{row.get('denominators', '')}` sum_minus_1=`{row.get('sum_minus_1', '')}`"
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge chunked anomaly CSV outputs into one resumable candidate dump.")
    parser.add_argument("root", type=Path, help="Directory containing chunk CSVs")
    parser.add_argument("--pattern", default="anomalies_N*_candidates_*.csv")
    parser.add_argument("--status-pattern", default="anomalies_N*_candidates_*.status.json")
    parser.add_argument("--out-csv", type=Path, default=Path("results/anomalies_candidates.csv"))
    parser.add_argument("--out-summary", type=Path, default=Path("results/anomaly_candidate_summary.md"))
    parser.add_argument("--N", type=int, default=None)
    parser.add_argument("--P", type=int, default=None)
    parser.add_argument("--low", default=None)
    parser.add_argument("--high", default=None)
    parser.add_argument("--expected-total", type=int, default=None)
    parser.add_argument("--expected-backbone-count", action="append", type=parse_expected_count, default=None)
    parser.add_argument("--expected-chunk-count", type=int, default=None)
    parser.add_argument("--denominator-max", type=int, default=65)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()

    paths = sorted(args.root.rglob(args.pattern))
    if not paths:
        raise SystemExit(f"No CSV files found under {args.root} matching {args.pattern}")
    status_paths = sorted(args.root.rglob(args.status_pattern))

    rows, fieldnames = load_rows(paths)
    statuses = load_statuses(status_paths)
    expected_counts = dict(args.expected_backbone_count or [])
    write_csv(args.out_csv, rows, fieldnames)
    write_summary(
        args.out_summary,
        rows,
        statuses,
        expected_counts,
        args.expected_total,
        args.expected_chunk_count,
        args.denominator_max,
        args.N,
        args.P,
        args.low,
        args.high,
    )
    print(f"Merged {len(paths)} files into {args.out_csv}")
    print(f"Status sidecars: {len(statuses)}")
    print(f"Rows: {len(rows)}")
    incomplete_chunks = [status for status in statuses if status.get("status") != "complete"]
    if incomplete_chunks:
        print("Incomplete chunks:")
        for status in incomplete_chunks:
            print(f"  {status.get('chunk_id', 'unknown')}: {status.get('status', 'unknown')}")
    if expected_counts:
        for backbone, expected_count in sorted(expected_counts.items()):
            actual = sum(1 for row in rows if backbone in row_backbones(row))
            status = "ok" if actual == expected_count else "incomplete"
            print(f"{backbone}: {actual}/{expected_count} {status}")
    else:
        backbone_counts = Counter(backbone for row in rows for backbone in row_backbones(row))
        for backbone, actual in sorted(backbone_counts.items()):
            print(f"{backbone}: {actual}")
    row_counts_ok = all(
        sum(1 for row in rows if backbone in row_backbones(row)) == expected_count
        for backbone, expected_count in expected_counts.items()
    )
    if args.require_complete:
        if not statuses:
            raise SystemExit("Completeness required, but no chunk status sidecars were found")
        if args.expected_chunk_count is not None and len(statuses) != args.expected_chunk_count:
            raise SystemExit(
                f"Completeness required, but found {len(statuses)} status sidecars instead of {args.expected_chunk_count}"
            )
        commit_shas = {str(status.get("commit_sha", "")) for status in statuses if status.get("commit_sha")}
        source_hashes = {str(status.get("cpp_sha256", "")) for status in statuses if status.get("cpp_sha256")}
        compile_flags = {str(status.get("compile_flags", "")) for status in statuses if status.get("compile_flags")}
        if len(commit_shas) > 1 or len(source_hashes) > 1 or len(compile_flags) > 1:
            raise SystemExit("Completeness required, but chunk build metadata is inconsistent")
        if incomplete_chunks:
            raise SystemExit("Completeness required, but one or more chunks are incomplete")
        if args.expected_total is not None and len(rows) != args.expected_total:
            raise SystemExit("Completeness required, but merged anomaly row count does not match expectation")
        if expected_counts and not row_counts_ok:
            raise SystemExit("Completeness required, but merged anomaly counts do not match expectations")


if __name__ == "__main__":
    main()
