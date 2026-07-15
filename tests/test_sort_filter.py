"""Tests for sort and filter commands."""

from __future__ import annotations

from tests.conftest import ShredtoolsRunner, load_mum_header


def test_sort_verify_unsorted_fails(runner: ShredtoolsRunner, bac_paths):
    result = runner.run("sort", str(bac_paths["bumbl"]), "--verify", "-s", "0", check=False)
    assert result.returncode == 1
    assert "not sorted" in result.stderr


def test_sort_verify_sorted_passes(runner: ShredtoolsRunner, bac_sorted):
    result = runner.run_ok("sort", str(bac_sorted), "--verify", "-s", "0")
    assert "sorted" in result.stderr


def test_sort_and_inplace_conflict(runner: ShredtoolsRunner, bac_paths, tmp_out):
    out = tmp_out / "sorted.bumbl"
    result = runner.run(
        "sort",
        str(bac_paths["bumbl"]),
        "-i",
        "-o",
        str(out),
        check=False,
    )
    assert result.returncode == 1
    assert "Cannot use both" in result.stderr


def test_sort_produces_sorted_output(runner: ShredtoolsRunner, bac_paths, tmp_out):
    out = tmp_out / "sorted.bumbl"
    runner.run_ok("sort", str(bac_paths["bumbl"]), "-s", "0", "-o", str(out))
    runner.run_ok("sort", str(out), "--verify", "-s", "0")
    header = load_mum_header(out)
    assert header["n_mums"] == 6848


def test_filter_reduces_or_preserves_mums(runner: ShredtoolsRunner, bac_sorted, tmp_out):
    out = tmp_out / "filtered.bumbl"
    runner.run_ok("filter", str(bac_sorted), "-o", str(out))
    header = load_mum_header(out)
    assert 0 < header["n_mums"] <= 6848
    assert header["blocks_stored"] is True
