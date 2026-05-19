#!/usr/bin/env python3
from __future__ import annotations

import csv
import itertools
import json
import re
import zipfile
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence


BACKBONE_DEFAULTS: tuple[tuple[str, str], ...] = (("2", "31"), ("19", "37"))


@dataclass(frozen=True)
class PrimeTopLayerProfile:
    prime: int
    max_valuation: int
    top_denominators: tuple[int, ...]
    residue_sum_mod_p: int

    @property
    def kills(self) -> bool:
        residue = self.residue_sum_mod_p - 1 if self.max_valuation == 0 else self.residue_sum_mod_p
        return residue % self.prime != 0


def normalize_header(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def extract_ints(text: str) -> list[int]:
    return [int(token) for token in re.findall(r"\d+", text)]


def parse_backbone_text(text: str) -> tuple[str, str]:
    numbers = extract_ints(text)
    if len(numbers) < 2:
        raise ValueError(f"Could not parse backbone from {text!r}")
    return tuple(str(value) for value in numbers[:2])  # type: ignore[return-value]


def parse_prime_labels(text: str) -> tuple[str, ...]:
    return tuple(str(value) for value in extract_ints(text))


def bits_to_labels(mask: int, labels: Sequence[str]) -> list[str]:
    return [label for index, label in enumerate(labels) if (mask >> index) & 1]


def labels_to_mask(selected: Iterable[str], labels: Sequence[str]) -> int:
    wanted = set(str(value) for value in selected)
    mask = 0
    for index, label in enumerate(labels):
        if label in wanted:
            mask |= 1 << index
    return mask


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


def find_uncovered_mask(mask_counts: dict[int, int], cover_mask: int) -> tuple[int, int] | None:
    for mask, count in sorted(mask_counts.items()):
        if mask & cover_mask == 0:
            return mask, count
    return None


def load_hyb_run(path: Path) -> dict:
    if path.suffix.lower() == ".zip":
        return load_hyb_run_from_zip(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    return extract_run(data, path)


def load_hyb_run_from_zip(path: Path) -> dict:
    with zipfile.ZipFile(path) as archive:
        json_members = [
            name
            for name in archive.namelist()
            if name.endswith("_hyb_cpp.json") and "progress" not in name and not name.endswith("_merged.json")
        ]
        if not json_members:
            raise ValueError(f"No *_hyb_cpp.json member found in {path}")
        if len(json_members) > 1:
            raise ValueError(f"Expected one hyb result JSON in {path}, found {json_members}")
        with archive.open(json_members[0]) as handle:
            data = json.load(handle)
    return extract_run(data, path)


def extract_run(payload: dict, origin: Path) -> dict:
    runs = payload.get("runs") or []
    if not runs:
        raise ValueError(f"{origin} has no runs[] payload")
    run = dict(runs[0])
    run["_completed"] = bool(payload.get("completed", False))
    run["_interrupted"] = bool(payload.get("interrupted", False))
    run["_origin"] = str(origin)
    return run


def get_mask_counts(run: dict) -> dict[int, int]:
    unique_masks = run.get("unique_masks") or {}
    return {int(mask): int(count) for mask, count in unique_masks.items()}


def get_test_labels(run: dict) -> list[str]:
    labels = run.get("tests")
    if not isinstance(labels, list) or not labels:
        raise ValueError("Run is missing tests[] labels")
    return [str(label) for label in labels]


def escape_summary(mask_counts: dict[int, int], labels: Sequence[str], backbone: tuple[str, str]) -> dict:
    cover_mask = labels_to_mask(backbone, labels)
    escaping_rows = []
    third_prime_counts: Counter[str] = Counter()
    overlap_masks: list[int] = []
    for mask, count in sorted(mask_counts.items()):
        if mask & cover_mask:
            continue
        kill_labels = bits_to_labels(mask, labels)
        escaping_rows.append({
            "mask": mask,
            "count": count,
            "kill_primes": kill_labels,
        })
        for prime in kill_labels:
            if prime not in backbone:
                third_prime_counts[prime] += count
        if not kill_labels:
            overlap_masks.append(mask)
    return {
        "backbone": backbone,
        "near_count": sum(row["count"] for row in escaping_rows),
        "unique_mask_count": len(escaping_rows),
        "escaping_rows": escaping_rows,
        "third_prime_counts": dict(sorted(third_prime_counts.items(), key=lambda item: (int(item[0]), item[1]))),
        "overlap_masks": overlap_masks,
    }


def pair_witness_rows(mask_counts: dict[int, int], labels: Sequence[str]) -> list[dict]:
    rows = []
    for left_index, right_index in itertools.combinations(range(len(labels)), 2):
        cover_mask = (1 << left_index) | (1 << right_index)
        witness = find_uncovered_mask(mask_counts, cover_mask)
        rows.append({
            "cover_labels": (labels[left_index], labels[right_index]),
            "cover_mask": cover_mask,
            "witness_mask": None if witness is None else witness[0],
            "witness_count": 0 if witness is None else witness[1],
            "witness_kill_primes": [] if witness is None else bits_to_labels(witness[0], labels),
            "passes": witness is None,
        })
    return rows


def mod_pow(base: int, exp: int, mod: int) -> int:
    result = 1
    value = base % mod
    power = exp
    while power > 0:
        if power & 1:
            result = (result * value) % mod
        value = (value * value) % mod
        power >>= 1
    return result


def mod_inverse_prime(value: int, prime: int) -> int:
    normalized = value % prime
    if normalized == 0:
        raise ValueError(f"{value} is not invertible modulo {prime}")
    return mod_pow(normalized, prime - 2, prime)


def prime_valuation(denominator: int, prime: int) -> int:
    value = denominator
    valuation = 0
    while value % prime == 0:
        value //= prime
        valuation += 1
    return valuation


def top_layer_profile(denominators: Sequence[int], prime: int) -> PrimeTopLayerProfile:
    if not denominators:
        raise ValueError("denominators must be non-empty")

    max_valuation = max(prime_valuation(denominator, prime) for denominator in denominators)
    top_denominators = tuple(denominator for denominator in denominators if prime_valuation(denominator, prime) == max_valuation)
    residue_sum = 0
    for denominator in top_denominators:
        reduced = denominator // (prime ** max_valuation)
        residue_sum = (residue_sum + mod_inverse_prime(reduced, prime)) % prime
    return PrimeTopLayerProfile(
        prime=prime,
        max_valuation=max_valuation,
        top_denominators=top_denominators,
        residue_sum_mod_p=residue_sum,
    )


def is_gap_leq_two(denominators: Sequence[int]) -> bool:
    if not denominators:
        return True
    return all((right - left) <= 2 for left, right in zip(denominators, denominators[1:]))


def gap_pattern(denominators: Sequence[int]) -> tuple[int, ...]:
    return tuple(right - left for left, right in zip(denominators, denominators[1:]))


def exact_sum_fraction(denominators: Sequence[int]) -> Fraction:
    total = Fraction(0, 1)
    for denominator in denominators:
        total += Fraction(1, denominator)
    return total


def format_fraction(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def parse_backbone_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"{path} has no CSV headers")
        header_map = {normalize_header(name): name for name in reader.fieldnames}

        backbone_col = choose_header(header_map, ["backbone", "escapedbackbone", "coverpair", "pair"])
        mask_col = choose_header(header_map, ["mask", "killmask", "maskdecimal", "maskint"])
        count_col = choose_header(header_map, ["count", "misscount", "nearcount", "escapingcount", "misses"], required=False)
        kills_col = choose_header(header_map, ["kills", "killprimes", "killlabels", "primes"], required=False)

        rows = []
        for row in reader:
            if not any((value or "").strip() for value in row.values()):
                continue
            backbone = parse_backbone_text(row[backbone_col])
            mask = int(extract_ints(row[mask_col])[0])
            count = int(extract_ints(row[count_col])[0]) if count_col and extract_ints(row[count_col]) else 1
            kill_primes = parse_prime_labels(row[kills_col]) if kills_col else ()
            rows.append({
                "backbone": backbone,
                "mask": mask,
                "count": count,
                "kill_primes": kill_primes,
            })
    return rows


def choose_header(header_map: dict[str, str], candidates: Sequence[str], required: bool = True) -> str | None:
    for candidate in candidates:
        if candidate in header_map:
            return header_map[candidate]
    if required:
        raise ValueError(f"Missing expected CSV headers from {sorted(header_map.values())}")
    return None


def compare_csv_to_summary(path: Path, summaries: Sequence[dict], labels: Sequence[str]) -> dict:
    rows = parse_backbone_csv(path)
    expected: dict[tuple[tuple[str, str], int], dict] = {}
    for summary in summaries:
        backbone = tuple(summary["backbone"])
        for row in summary["escaping_rows"]:
            expected[(backbone, row["mask"])] = {
                "count": row["count"],
                "kill_primes": tuple(row["kill_primes"]),
            }

    seen_keys = set()
    mismatches = []
    for row in rows:
        key = (tuple(row["backbone"]), row["mask"])
        seen_keys.add(key)
        expected_row = expected.get(key)
        if expected_row is None:
            mismatches.append(f"Unexpected CSV row for backbone={row['backbone']} mask={row['mask']}")
            continue
        if row["count"] != expected_row["count"]:
            mismatches.append(
                f"Count mismatch for backbone={row['backbone']} mask={row['mask']}: "
                f"csv={row['count']} json={expected_row['count']}"
            )
        if row["kill_primes"] and tuple(row["kill_primes"]) != expected_row["kill_primes"]:
            mismatches.append(
                f"Kill-prime mismatch for backbone={row['backbone']} mask={row['mask']}: "
                f"csv={row['kill_primes']} json={expected_row['kill_primes']}"
            )

    missing = sorted(set(expected) - seen_keys)
    for backbone, mask in missing:
        mismatches.append(f"Missing CSV row for backbone={backbone} mask={mask}")

    return {
        "row_count": len(rows),
        "expected_row_count": len(expected),
        "mismatches": mismatches,
    }
