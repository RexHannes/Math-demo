# Erdős #287: modulus-33 regime probe

This package investigates the composite-modulus regime around

```text
33 = 3 * 11
```

for the gap-2 Egyptian-fraction interval certificate pipeline.

The key local observation from preliminary runs is:

```text
Modulo 33 alone blocks every feasible support interval [m,M]
for 242 <= M <= 362.
```

Here "feasible" means `2 <= m <= m_star(M)`, where `m_star(M)` is the smallest `m` such that the full harmonic interval sum `sum_{n=m}^M 1/n` can still reach 1. If `m > m_star(M)`, the full interval already has reciprocal mass below 1, so no subset can work.

This is stronger than the earlier pipeline statement that modulus 33 starts appearing around `M ≈ 245`: the pipeline only selects 33 after direct prime or linked-pair certificates fail. The pure modulus-33 test shows the 33 obstruction itself is already active from `M=242` through `M=362`.

## Local run

```bash
g++ -O3 -std=c++17 src/probe_specific_modulus.cpp -o probe_modulus
mkdir -p results
./probe_modulus \
  --M-min 200 \
  --M-max 450 \
  --modulus 33 \
  --prefix results/mod33_M200_450
```

Outputs:

```text
results/mod33_M200_450_summary_by_M.csv
results/mod33_M200_450_reachable_records.csv
results/mod33_M200_450_all_records.csv
```

## GitHub Actions

Two workflows are included:

```text
.github/workflows/erdos_mod33_smalltest.yml
.github/workflows/erdos_mod33_regime.yml
```

Run the small test first, then the larger matrix scan.

## Research goal

The current forced-pair route remains parity-linked because it relies on adjacent prime pairs such as `(106,107)`. The modulus-33 regime is interesting because it is a composite-modulus obstruction that does not obviously require new Sophie-Germain or safe-prime input.

The next questions are:

1. Can the `33` regime be proved cleanly by hand for `242 <= M <= 362`?
2. What replaces `33` after `M=363`?
3. Are there later composite regimes similar to `33`, possibly involving higher powers or layers of `3` and `11`?
4. Does a broader family of composite moduli progressively absorb intervals that otherwise need linked-prime forced-pair certificates?
