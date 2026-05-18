#!/usr/bin/env python3
"""Merge hyb cover chunk JSON files and compute the combined hitting set."""
from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path
from typing import Iterable


def bits_to_labels(mask: int, labels: list[str]) -> list[str]:
    return [label for index, label in enumerate(labels) if (mask >> index) & 1]


def minimal_hitting_sets(unique_masks: Iterable[int], num_tests: int, max_solutions: int = 20) -> tuple[int, list[int]]:
    masks = sorted(set(unique_masks), key=lambda value: (value.bit_count(), value))
    if not masks:
        return 0, [0]

    solutions: list[int] = []

    def search(start_bit: int, remaining: int, cover_mask: int) -> bool:
        if remaining == 0:
            if all(candidate_mask & cover_mask for candidate_mask in masks):
                solutions.append(cover_mask)
                return len(solutions) >= max_solutions
            return False
        for bit in range(start_bit, num_tests - remaining + 1):
            if search(bit + 1, remaining - 1, cover_mask | (1 << bit)):
                return True
        return False

    for size in range(1, num_tests + 1):
        solutions.clear()
        search(0, size, 0)
        if solutions:
            return size, solutions
    return num_tests + 1, []


def load_run(path: Path) -> dict:
    data = json.loads(path.read_text())
    runs = data.get("runs") or []
    if not runs:
        raise ValueError(f"{path} has no runs")
    return runs[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path, help="Directory containing downloaded/extracted hyb JSON artifacts")
    parser.add_argument("--out", type=Path, default=Path("results/merged_hyb_cover.json"))
    args = parser.parse_args()

    files = sorted(args.root.rglob("*_hyb_cpp.json"))
    if not files:
        raise SystemExit(f"No *_hyb_cpp.json files found under {args.root}")

    mask_counts: Counter[int] = Counter()
    chunk_summaries = []
    tests: list[str] | None = None
    completed = True
    total_nodes = 0
    total_near = 0
    first_uncovered = None

    for path in files:
        run = load_run(path)
        if tests is None:
            tests = list(run["tests"])
        elif tests != run["tests"]:
            raise ValueError(f"Test labels differ in {path}")

        unique_masks = run.get("unique_masks") or {}
        mask_counts.update({int(mask): int(count) for mask, count in unique_masks.items()})
        total_nodes += int(run.get("nodes", 0))
        total_near += int(run.get("near_count", 0))
        completed = completed and bool(json.loads(path.read_text()).get("completed", False))
        if first_uncovered is None and run.get("first_uncovered") is not None:
            first_uncovered = run["first_uncovered"]
        chunk_summaries.append({
            "file": str(path),
            "completed": bool(json.loads(path.read_text()).get("completed", False)),
            "start_min": run.get("start_min"),
            "start_max": run.get("start_max"),
            "nodes": run.get("nodes", 0),
            "near_count": run.get("near_count", 0),
            "unique_mask_count": len(unique_masks),
        })

    assert tests is not None
    min_size, cover_masks = minimal_hitting_sets(mask_counts.keys(), len(tests))
    top_masks = sorted(mask_counts.items(), key=lambda item: (-item[1], item[0]))[:20]

    output = {
        "generated_at": str(int(time.time())),
        "completed": completed,
        "chunk_count": len(chunk_summaries),
        "nodes": total_nodes,
        "near_count": total_near,
        "unique_mask_count": len(mask_counts),
        "minimal_cover_size": min_size if completed else None,
        "partial_minimal_cover_size_so_far": min_size,
        "minimal_covers_sample": [bits_to_labels(mask, tests) for mask in cover_masks],
        "first_uncovered": first_uncovered,
        "top_masks": [
            {"mask": mask, "count": count, "kills": bits_to_labels(mask, tests)}
            for mask, count in top_masks
        ],
        "chunks": chunk_summaries,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
