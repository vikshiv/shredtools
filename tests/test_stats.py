"""Tests for shredtools stats."""

from __future__ import annotations

import json

from tests.conftest import ShredtoolsRunner, parse_stats_json


def test_stats_bac_human(runner: ShredtoolsRunner, bac_paths):
    result = runner.run_ok("stats", str(bac_paths["bumbl"]))
    assert "num_mums: 6848" in result.stdout
    assert "num_seqs: 5" in result.stdout
    assert "lengths:" in result.stdout


def test_stats_json(runner: ShredtoolsRunner, bac_paths):
    report = parse_stats_json(runner, bac_paths["bumbl"])
    assert report["format"] == "bumbl"
    assert report["num_mums"] == 6848
    assert report["num_seqs"] == 5
    assert "lengths" in report["associated_files"]


def test_stats_missing_file(runner: ShredtoolsRunner):
    result = runner.run("stats", "/no/such/file.bumbl", check=False)
    assert result.returncode == 1
    assert "not found" in result.stderr
