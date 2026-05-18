# Hyb Cover Speedups

## Current safe path

For `N=70`, prefer the chunked C++ workflow:

```text
Erdős 287 hyb cover chunks (C++)
```

It runs these starting-denominator chunks in parallel:

```text
2-10
11-20
21-30
31-40
41-50
51-70
```

Each chunk writes the full `unique_masks` dictionary into its `hyb_cpp.json` artifact. The merge job combines those dictionaries and computes the hitting set from the merged masks.

Use:

```text
N=70
P=70
low=0.999
high=1.001
shell_timeout_minutes=75
progress_every=10000000
```

If every chunk finishes, the merged artifact can certify `minimal_cover_size`. If one or more chunks time out, the merged artifact still reports `partial_minimal_cover_size_so_far`, but the final `minimal_cover_size` remains `null`.

## Meet in the middle note

Meet in the middle is promising, but it must not combine partial kill masks with:

```cpp
combined_mask = left_mask | right_mask;
```

That shortcut is not mathematically valid here because modular residue defects can cancel across the two halves.

A correct meet-in-the-middle implementation must store enough per-prime p-adic state for each half, including the maximum valuation and top-layer residue contribution. At join time it must combine those states and recompute the hybrid mask for the combined set.

Safe implementation rule:

```cpp
combined_state = combine(left_state, right_state);
combined_mask = compute_mask(combined_state);
```

Before trusting any meet-in-the-middle rewrite, compare it against the existing DFS for small values such as `N=25`, `N=30`, and `N=35`. The unique masks and minimal cover size must match exactly.
