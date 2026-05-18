# Erdős composite-modulus regime scan

This track runs a broader GitHub Actions scan for composite moduli that may produce finite full-block windows, similar to the known mod-33 window.

It does not use CP-SAT. The method is still the direct fixed-modulus no-adjacent-skip residue automaton:

```text
choose modulus R
-> scan an M-range
-> detect full-block windows
-> inspect mechanism only if something appears
```

## Workflows

- `.github/workflows/erdos_multi_moduli_regime_scan.yml`
- `.github/workflows/erdos_single_modulus_examples.yml`

## Recommended first run

Run `Erdős multi-moduli regime scan` with defaults:

```text
M_min = 200
M_max = 700
moduli = 33 39 55 77 99 121 143 165 187 231 297 363
keep_examples = false
```

This is the broad scan. It keeps outputs smaller by omitting example residue lists during the first pass.

## After the run

Download all workflow artifacts into one folder, then run:

```bash
python3 scripts/summarize_multi_moduli.py artifacts --min-len 20
```

The main outputs are:

- `multi_moduli_summary_by_modulus.csv`
- `multi_moduli_summary_windows_minlen20.csv`
- `multi_moduli_summary_all_by_M.csv`

The most important file is:

- `multi_moduli_summary_windows_minlen20.csv`

That tells us whether any modulus has a long full-block window.

## Follow-up

If the broad scan finds a promising window for modulus `R`, run:

- `Erdős single modulus with examples`

using:

```text
M_min = A
M_max = B
modulus = R
```

where `[A, B]` is the window to inspect. That keeps the nonzero residue weights so the mechanism can be examined directly.

## Interpretation

- One isolated window: useful finite example, but probably a local episode.
- Several long windows with overlapping coverage: more serious composite-regime evidence.
- A visible family of windows with coherent structure: worth much deeper follow-up.
