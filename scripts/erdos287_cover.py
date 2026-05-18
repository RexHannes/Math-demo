#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
import time
from typing import Callable, Dict, List, Sequence, Tuple


def primes_upto(n: int) -> List[int]:
    if n < 2:
        return []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            start = i * i
            step = i
            sieve[start:n + 1:step] = [False] * (((n - start) // step) + 1)
    return [i for i, ok in enumerate(sieve) if ok]


def prime_powers_upto(n: int) -> List[Tuple[int, int, int]]:
    out = []
    for p in primes_upto(n):
        q = p
        r = 1
        while q <= n:
            out.append((p, r, q))
            r += 1
            q *= p
    return out


def suffix_harmonic_float(N: int) -> List[float]:
    suffix = [0.0] * (N + 3)
    for k in range(N, 1, -1):
        suffix[k] = suffix[k + 1] + 1.0 / k
    return suffix


def lcm_many(nums: Sequence[int]) -> int:
    L = 1
    for n in nums:
        L = L * n // math.gcd(L, n)
    return L


def cleared_numerator(S: Sequence[int]) -> int:
    L = lcm_many(S)
    return sum(L // n for n in S) - L


def mask_for_A(A: int, tests: Sequence[int]) -> int:
    mask = 0
    for i, q in enumerate(tests):
        if A % q != 0:
            mask |= 1 << i
    return mask


def minimal_hitting_sets(unique_masks: Sequence[int], num_tests: int, max_solutions: int = 20) -> Tuple[int, List[int]]:
    if not unique_masks:
        return 0, [0]

    masks = sorted(set(unique_masks), key=lambda x: (x.bit_count(), x))
    all_bits = list(range(num_tests))

    for k in range(1, num_tests + 1):
        solutions = []
        for combo_bits in itertools.combinations(all_bits, k):
            cover_mask = 0
            for bit in combo_bits:
                cover_mask |= 1 << bit
            if all((candidate_mask & cover_mask) != 0 for candidate_mask in masks):
                solutions.append(cover_mask)
                if len(solutions) >= max_solutions:
                    return k, solutions
        if solutions:
            return k, solutions

    return num_tests + 1, []


def bits_to_labels(mask: int, labels: Sequence[str]) -> List[str]:
    return [labels[i] for i in range(len(labels)) if (mask >> i) & 1]


def write_progress_snapshot(progress_out: str | None, payload: Dict) -> None:
    if not progress_out:
        return
    path = Path(progress_out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def enumerate_gap2_near(
    N: int,
    low: float,
    high: float,
    test_moduli: Sequence[int],
    sample_limit: int = 5,
    progress_every: int = 0,
    progress_callback: Callable[[Dict], None] | None = None,
) -> Dict:
    suffix = suffix_harmonic_float(N)
    nodes = 0
    near_count = 0
    exact_hits: List[List[int]] = []
    sample_candidates = []
    unique_masks: Dict[int, int] = {}
    first_uncovered = None
    started = time.time()

    def visit(S: List[int], last: int, current_sum: float) -> None:
        nonlocal nodes, near_count, first_uncovered
        nodes += 1

        if progress_every and nodes % progress_every == 0:
            payload = {
                "N": N,
                "nodes": nodes,
                "near_count": near_count,
                "last": last,
                "sum_float": current_sum,
                "elapsed_seconds": time.time() - started,
                "first_uncovered": first_uncovered,
            }
            print(
                f"progress N={N}: nodes={nodes:,}, near={near_count:,}, "
                f"last={last}, sum={current_sum:.8f}",
                flush=True,
            )
            if progress_callback:
                progress_callback(payload)

        if current_sum > high:
            return

        if last + 1 <= N and current_sum + suffix[last + 1] < low:
            return

        if low <= current_sum <= high:
            A = cleared_numerator(S)
            if A == 0:
                exact_hits.append(list(S))

            kill_mask = mask_for_A(A, test_moduli)
            unique_masks[kill_mask] = unique_masks.get(kill_mask, 0) + 1
            near_count += 1

            if kill_mask == 0 and first_uncovered is None:
                first_uncovered = list(S)

            if len(sample_candidates) < sample_limit:
                sample_candidates.append({
                    "S": list(S),
                    "sum_float": current_sum,
                    "A": A,
                    "kill_mask": kill_mask
                })

        for nxt in (last + 1, last + 2):
            if nxt <= N:
                S.append(nxt)
                visit(S, nxt, current_sum + 1.0 / nxt)
                S.pop()

    for start in range(2, N + 1):
        visit([start], start, 1.0 / start)

    return {
        "N": N,
        "low": low,
        "high": high,
        "nodes": nodes,
        "near_count": near_count,
        "exact_hits": exact_hits[:sample_limit],
        "exact_hit_count": len(exact_hits),
        "unique_mask_count": len(unique_masks),
        "unique_masks": unique_masks,
        "first_uncovered": first_uncovered,
        "sample_candidates": sample_candidates,
        "seconds": time.time() - started
    }


def run_one(N: int, low: float, high: float, P: int, mode: str) -> Dict:
    if mode == "hyb":
        test_moduli = primes_upto(P)
        labels = [str(p) for p in test_moduli]
    elif mode == "lift":
        prime_powers = prime_powers_upto(P)
        test_moduli = [q for _, _, q in prime_powers]
        labels = [f"{p}^{r}" if r > 1 else str(p) for p, r, _ in prime_powers]
    else:
        raise ValueError("mode must be hyb or lift")

    result = enumerate_gap2_near(
        N,
        low,
        high,
        test_moduli,
        progress_every=run_one.progress_every,
        progress_callback=run_one.progress_callback,
    )
    unique_masks = list(result["unique_masks"].keys())
    min_size, cover_masks = minimal_hitting_sets(unique_masks, len(test_moduli))
    cover_labels = [bits_to_labels(mask, labels) for mask in cover_masks]

    top_masks = sorted(
        result["unique_masks"].items(),
        key=lambda kv: kv[1],
        reverse=True
    )[:20]

    return {
        "mode": mode,
        "N": N,
        "P": P,
        "low": low,
        "high": high,
        "test_count": len(test_moduli),
        "tests": labels,
        "nodes": result["nodes"],
        "near_count": result["near_count"],
        "exact_hit_count": result["exact_hit_count"],
        "exact_hits_sample": result["exact_hits"],
        "unique_mask_count": result["unique_mask_count"],
        "minimal_cover_size": min_size,
        "minimal_covers_sample": cover_labels,
        "first_uncovered": result["first_uncovered"],
        "top_masks": [
            {
                "count": count,
                "kills": bits_to_labels(mask, labels),
                "mask": mask
            }
            for mask, count in top_masks
        ],
        "sample_candidates": [
            {
                **candidate,
                "kills": bits_to_labels(candidate["kill_mask"], labels)
            }
            for candidate in result["sample_candidates"]
        ],
        "seconds": result["seconds"]
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, required=True)
    parser.add_argument("--P", type=int, default=None, help="max prime / prime-power bound; default N")
    parser.add_argument("--low", type=float, default=0.99)
    parser.add_argument("--high", type=float, default=1.01)
    parser.add_argument("--mode", choices=["hyb", "lift", "both"], default="both")
    parser.add_argument("--out", type=str, default=None)
    parser.add_argument("--progress-every", type=int, default=100_000)
    parser.add_argument("--progress-out", type=str, default=None)
    args = parser.parse_args()

    P = args.P if args.P is not None else args.N
    modes = ["hyb", "lift"] if args.mode == "both" else [args.mode]
    outputs = []

    for mode in modes:
        mode_progress_out = None
        if args.progress_out:
            base = Path(args.progress_out)
            suffix = base.suffix or ".json"
            if args.mode == "both":
                mode_progress_out = str(base.with_name(f"{base.stem}_{mode}{suffix}"))
            else:
                mode_progress_out = str(base)

        run_one.progress_every = args.progress_every
        run_one.progress_callback = lambda payload, mode=mode, P=P, low=args.low, high=args.high, out=mode_progress_out: write_progress_snapshot(out, {
            "mode": mode,
            "P": P,
            "low": low,
            "high": high,
            **payload,
        })

        print(f"Running mode={mode}, N={args.N}, P={P}, window=[{args.low},{args.high}]")
        result = run_one(args.N, args.low, args.high, P, mode)
        outputs.append(result)
        print(json.dumps({
            "mode": result["mode"],
            "N": result["N"],
            "P": result["P"],
            "nodes": result["nodes"],
            "near_count": result["near_count"],
            "exact_hit_count": result["exact_hit_count"],
            "unique_mask_count": result["unique_mask_count"],
            "minimal_cover_size": result["minimal_cover_size"],
            "minimal_covers_sample": result["minimal_covers_sample"][:5],
            "first_uncovered": result["first_uncovered"],
            "seconds": result["seconds"]
        }, indent=2))

    final = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "runs": outputs
    }

    if args.progress_out:
        write_progress_snapshot(args.progress_out, {
            "generated_at": final["generated_at"],
            "completed": True,
            "runs": [
                {
                    "mode": run["mode"],
                    "N": run["N"],
                    "P": run["P"],
                    "nodes": run["nodes"],
                    "near_count": run["near_count"],
                    "unique_mask_count": run["unique_mask_count"],
                    "minimal_cover_size": run["minimal_cover_size"],
                    "first_uncovered": run["first_uncovered"],
                    "seconds": run["seconds"],
                }
                for run in outputs
            ],
        })

    if args.out:
        with open(args.out, "w", encoding="utf-8") as handle:
            json.dump(final, handle, indent=2)
        print(f"Wrote {args.out}")
    else:
        print(json.dumps(final, indent=2))


if __name__ == "__main__":
    main()
