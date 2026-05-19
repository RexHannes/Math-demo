# Input Artifacts

Place uploaded artifact files here so Codex or GitHub Actions can read them from disk instead of from chat.

Recommended bundle for the N65 anomaly task:

- `erdos287-cover-hyb-cpp-N65-P65-S2-0.zip`
  Main file. Expected to contain the full N65 `unique_masks` dictionary and cover summary.
- `N65_backbone_escaping_masks.csv`
  Helper file. Derived mask-level summary for the `{2,31}` and `{19,37}` backbone checks.
- `erdos287-cover-hyb-cpp-N70-P70-merged.zip`
- `erdos287-cover-hyb-cpp-N70-P70-merged (1).zip`
  Background-only merged summaries.

Expected workflow:

1. Inspect the N65 ZIP and CSV first.
2. Write scripts that parse files from `data/input/`.
3. Do not print or paste large artifact contents.
4. If candidate-level anomaly inspection is needed, use the C++ source in this repo and add a selective dump mode instead of relying on the ZIP alone.
