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
