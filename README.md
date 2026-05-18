# Erdős interval diagnostic: GitHub Actions batch runner

This repository scaffold runs the high-modulus interval-certificate diagnostic for the Erdős #287-style gap≤2 Egyptian-fraction problem.

Main experiment:

```text
mixed_limit = 10000
M = 101..200
dense_min_len = 20
```

## Files

- `src/analyze_certificate_dependencies.cpp` — C++ diagnostic program.
- `.github/workflows/erdos_mixed10000_smalltest.yml` — small manual test workflow.
- `.github/workflows/erdos_mixed10000.yml` — matrix workflow for `M=101..200` in 10 chunks.
- `scripts/summarize_artifacts.py` — local script to merge downloaded artifacts.

## How to run

1. Create a public GitHub repository.
2. Upload these files, preserving paths.
3. Go to **Actions**.
4. First run **Erdős mixed-modulus small test** with:

```text
start_M = 101
end_M = 102
mixed_limit = 10000
prime_limit = 1000
dense_min_len = 20
J = 12
```

5. If that passes, run **Erdős mixed-modulus diagnostic** with defaults.

The matrix workflow splits the calculation into `101-110`, `111-120`, ..., `191-200`, with `max-parallel: 4`. Each job runs one `M` at a time and writes separate `results/diag_M...` files, so partial progress is kept even if a job times out.

## After the run

Download all artifacts from the workflow run, unzip them into a directory named `artifacts`, then run:

```bash
python scripts/summarize_artifacts.py artifacts
```

The decisive columns are:

```text
mixed_or_large_modulus
dense_forced_holes
linked_prime_pair
survivor
```

The main research question is whether raising `mixed_limit` to `10000` makes the `linked_prime_pair` / forced-hole dependency collapse for `M <= 200`.

## Notes

The program uses `--prefix`, not `--out`, and writes:

```text
<prefix>_records.csv
<prefix>_summary_by_M.csv
<prefix>_report.json
```

## Mod-33 Probe

The repo now also contains a separate composite-modulus probe centered on the `33 = 3 * 11` regime:

- `src/probe_specific_modulus.cpp`
- `.github/workflows/erdos_mod33_smalltest.yml`
- `.github/workflows/erdos_mod33_regime.yml`
- `scripts/summarize_modulus_probe.py`
- `preliminary_results/probe33_230_370_*.csv`
- `docs/mod33_regime_probe.md`

This line is worth trying. It is not a replacement for the original mixed-modulus pipeline, but it pushes in a genuinely different direction: whether a clean composite modulus can absorb part of the range without leaning on linked-prime adjacent-pair certificates.

Suggested order:

1. Run `Erdős mod-33 small test`.
2. If it passes, run `Erdős mod-33 regime scan`.
3. Download artifacts and summarize them with:

```bash
python3 scripts/summarize_modulus_probe.py artifacts
```

## Multi-Moduli Scan

The repo also now supports the broader GitHub Actions scan discussed for testing whether mod `33` is isolated or part of a larger composite-modulus pattern.

- `.github/workflows/erdos_multi_moduli_regime_scan.yml`
- `.github/workflows/erdos_single_modulus_examples.yml`
- `scripts/summarize_multi_moduli.py`
- `docs/multi_moduli_regime_scan.md`

Recommended first run:

```text
M_min = 200
M_max = 700
moduli = 33 39 55 77 99 121 143 165 187 231 297 363
keep_examples = false
```

Then summarize downloaded artifacts with:

```bash
python3 scripts/summarize_multi_moduli.py artifacts --min-len 20
```

Only if a modulus shows a promising full-block window should you run the single-modulus examples workflow for mechanism inspection.

## Hybrid / Lift Covers

The repo also now includes a separate Python-based experiment for the stronger certificate hierarchy:

- `scripts/erdos287_cover.py`
- `.github/workflows/erdos287-cover.yml`

This computes the near-candidate cover statistics for:

- `C_hyb(N)`: prime-modulus hybrid cover
- `C_lift(N)`: prime-power lifting cover

Recommended order:

```bash
python3 scripts/erdos287_cover.py --N 50 --P 50 --mode both --out results/N50_P50.json
python3 scripts/erdos287_cover.py --N 60 --P 60 --mode both --out results/N60_P60.json
python3 scripts/erdos287_cover.py --N 70 --P 70 --mode both --out results/N70_P70.json
```

For GitHub Actions, run `Erdős 287 cover tests` with modest values first, then scale up to `N=90` only after the smaller jobs finish cleanly.

Operational note: the Python workflow now emits progress heartbeats, writes progress snapshots into `results/`, defaults to `hyb`, and uses a shell timeout so logs and partial outputs can still upload even when the full computation is too slow.

## Faster C++ Hybrid Cover

For heavier `hyb` runs, the repo also includes:

- `src/erdos287_cover_hyb.cpp`
- `.github/workflows/erdos287-cover-hyb-cpp.yml`

This C++ path avoids the Python brute-force overhead and is the better next move for jobs like:

- `N = 55, P = 55, low = 0.99, high = 1.01`
- `N = 60, P = 60, low = 0.999, high = 1.001`

It still uses the same gap-2 enumeration idea, so it is not a magic cure for exponential growth, but it should be much more practical for the next round of `hyb` experiments.
