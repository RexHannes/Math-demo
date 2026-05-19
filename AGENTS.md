# AGENTS.md

This repo contains Erdős 287 interval / cover experiments, including the hybrid C++ cover scanner.

General rules:
- Do not paste large JSON, CSV, log, or ZIP contents into chat.
- Read uploaded artifacts from disk, preferably from `data/input/`.
- Keep claims conservative: distinguish mask-level verification from candidate-level structural analysis.

For the N65 anomaly task:
- Main artifact: `data/input/erdos287-cover-hyb-cpp-N65-P65-S2-0.zip`
- Helper artifact: `data/input/N65_backbone_escaping_masks.csv`
- Background only: the two `N70` merged ZIPs

Important source paths already in this repo:
- `src/erdos287_cover_hyb.cpp`
- `scripts/erdos287_cover.cpp`

Interpretation rules:
- If only ZIP/CSV artifacts are available, only mask-level verification is possible.
- Candidate-level anomaly analysis requires rerunning the search with source code and a selective dump mode.

Preferred outputs:
- Put reproducible scripts in `scripts/`.
- Put generated summaries or CSVs in `results/`.
- Avoid noisy one-off notebooks unless explicitly requested.

Verification:
- If Python scripts change, run the relevant Python command or tests.
- If C++ scanner code changes, compile it before finalizing.
