from __future__ import annotations

import json
import sys
import zipfile
from pathlib import Path
import subprocess


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import mask_analysis_common as mac


FIXTURE_PAYLOAD = {
    "generated_at": "1779134030",
    "completed": True,
    "interrupted": False,
    "runs": [
        {
            "mode": "hyb",
            "N": 20,
            "P": 20,
            "low": 0.999,
            "high": 1.001,
            "test_count": 8,
            "tests": ["2", "3", "5", "7", "11", "13", "17", "19"],
            "nodes": 2296,
            "near_count": 21,
            "exact_hit_count": None,
            "exact_hits_sample": [],
            "unique_mask_count": 9,
            "unique_masks": {"127": 1, "191": 2, "245": 1, "247": 1, "249": 1, "251": 4, "253": 2, "254": 1, "255": 8},
            "minimal_cover_size": 1,
            "minimal_covers_sample": [["11"], ["13"]],
            "first_uncovered": None,
            "top_masks": [],
            "sample_candidates": [],
            "seconds": 0.0007603340,
        }
    ],
}


def write_fixture(path: Path) -> Path:
    path.write_text(json.dumps(FIXTURE_PAYLOAD), encoding="utf-8")
    return path


def test_load_hyb_run_from_json_fixture(tmp_path: Path) -> None:
    run = mac.load_hyb_run(write_fixture(tmp_path / "fixture.json"))
    assert run["near_count"] == 21
    assert run["minimal_cover_size"] == 1
    assert len(mac.get_mask_counts(run)) == 9


def test_load_hyb_run_from_zip(tmp_path: Path) -> None:
    zip_path = tmp_path / "fixture.zip"
    payload = json.dumps(FIXTURE_PAYLOAD)
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr("results/N55_P55_S2_0_hyb_cpp.json", payload)
        archive.writestr("results/N55_P55_S2_0_hyb_cpp_progress.json", "{}")
    run = mac.load_hyb_run(zip_path)
    assert run["P"] == 20
    assert run["N"] == 20


def test_minimal_hitting_set_reproduces_fixture(tmp_path: Path) -> None:
    run = mac.load_hyb_run(write_fixture(tmp_path / "fixture.json"))
    labels = mac.get_test_labels(run)
    mask_counts = mac.get_mask_counts(run)
    min_size, solutions = mac.minimal_hitting_sets(mask_counts, len(labels))
    assert min_size == 1
    assert any(set(mac.bits_to_labels(mask, labels)) == {"11"} for mask in solutions)


def test_escape_summary_for_fixture(tmp_path: Path) -> None:
    run = mac.load_hyb_run(write_fixture(tmp_path / "fixture.json"))
    labels = mac.get_test_labels(run)
    mask_counts = mac.get_mask_counts(run)
    summary = mac.escape_summary(mask_counts, labels, ("3", "7"))
    assert summary["near_count"] > 0
    assert summary["unique_mask_count"] > 0
    for row in summary["escaping_rows"]:
        kills = set(row["kill_primes"])
        assert "3" not in kills
        assert "7" not in kills


def test_helper_csv_crosscheck(tmp_path: Path) -> None:
    run = mac.load_hyb_run(write_fixture(tmp_path / "fixture.json"))
    labels = mac.get_test_labels(run)
    mask_counts = mac.get_mask_counts(run)
    summaries = [mac.escape_summary(mask_counts, labels, ("3", "7"))]
    csv_path = tmp_path / "escaping.csv"
    csv_path.write_text(
        "backbone,mask,count,kill_primes\n"
        f"\"{{3,7}}\",{summaries[0]['escaping_rows'][0]['mask']},{summaries[0]['escaping_rows'][0]['count']},"
        f"\"{';'.join(summaries[0]['escaping_rows'][0]['kill_primes'])}\"\n",
        encoding="utf-8",
    )
    comparison = mac.compare_csv_to_summary(csv_path, summaries, labels)
    assert comparison["mismatches"] == []


def test_top_layer_profile_matches_expected_behavior() -> None:
    profile = mac.top_layer_profile([4, 6, 8, 9], 2)
    assert profile.max_valuation == 3
    assert profile.top_denominators == (8,)
    assert profile.residue_sum_mod_p == 1
    assert profile.kills is True


def test_gap_helpers() -> None:
    assert mac.is_gap_leq_two([5, 6, 8, 9, 11]) is True
    assert mac.is_gap_leq_two([5, 8]) is False
    assert mac.gap_pattern([5, 6, 8, 9, 11]) == (1, 2, 1, 2)


def test_exact_sum_fraction() -> None:
    total = mac.exact_sum_fraction([2, 3, 6])
    assert total == 1
    assert mac.format_fraction(total) == "1/1"


def test_merge_anomaly_chunk_csvs(tmp_path: Path) -> None:
    first = tmp_path / "anomalies_N65_candidates_S2_10.csv"
    second = tmp_path / "anomalies_N65_candidates_S11_20.csv"
    header = "candidate_id,escaped_backbone,denominators,length,kill_primes,contains_37_38,contains_61_62,contains_adjacent_factor_pair,sum_minus_1\n"
    first.write_text(
        header + "1,2+31,5;6;7,3,3;5,false,false,false,1/100\n",
        encoding="utf-8",
    )
    second.write_text(
        header + "2,19+37,7;8;9,3,2;13,true,false,true,1/200\n",
        encoding="utf-8",
    )

    out_csv = tmp_path / "merged.csv"
    out_summary = tmp_path / "summary.md"
    subprocess.run(
        [
            "python3",
            str(Path(__file__).resolve().parents[1] / "scripts" / "merge_anomaly_chunk_csvs.py"),
            str(tmp_path),
            "--out-csv",
            str(out_csv),
            "--out-summary",
            str(out_summary),
        ],
        check=True,
    )

    merged = out_csv.read_text(encoding="utf-8")
    assert "S2_10:1" in merged
    assert "S11_20:2" in merged
    summary = out_summary.read_text(encoding="utf-8")
    assert "anomaly_rows: `2`" in summary
    assert "`2+31`: `1`" in summary
    assert "`19+37`: `1`" in summary
