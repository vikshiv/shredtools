"""Tests for shredtools extract."""

from __future__ import annotations

import pytest

from shredtools import utils as sutils
from shredtools.extract_from_mums import get_mum_ranges_flanks
from tests.conftest import ShredtoolsRunner

BAC_REGION = "GCF_000009045.1_ASM904v1_genomic.fna:100000-200000"


def test_extract_bed(
    runner: ShredtoolsRunner,
    bac_sorted,
    bac_multi_index,
    bac_paths,
    tmp_out,
):
    prefix = tmp_out / "extract"
    result = runner.run_ok(
        "extract",
        str(bac_sorted),
        "-s",
        "0",
        "-r",
        BAC_REGION,
        "-b",
        str(bac_multi_index),
        "-l",
        str(bac_paths["lengths"]),
        "-o",
        str(prefix),
    )
    bed = tmp_out / "extract.bed"
    assert bed.is_file()
    lines = bed.read_text().strip().splitlines()
    assert len(lines) == 5  # one row per sequence
    for line in lines:
        parts = line.split("\t")
        assert len(parts) == 4
        assert int(parts[1]) < int(parts[2])
    assert "left margin" in result.stderr
    assert "right margin" in result.stderr


def test_extract_stdout_bed(
    runner: ShredtoolsRunner,
    bac_sorted,
    bac_multi_index,
    bac_paths,
):
    result = runner.run_ok(
        "extract",
        str(bac_sorted),
        "-s",
        "0",
        "-r",
        BAC_REGION,
        "-b",
        str(bac_multi_index),
        "-l",
        str(bac_paths["lengths"]),
    )
    assert result.stdout.count("\t") >= 4


def test_extract_sequence_subset(
    runner: ShredtoolsRunner,
    bac_sorted,
    bac_multi_index,
    bac_paths,
    tmp_out,
):
    prefix = tmp_out / "extract_x"
    runner.run_ok(
        "extract",
        str(bac_sorted),
        "-s",
        "0",
        "-r",
        BAC_REGION,
        "-x",
        "0",
        "2",
        "-b",
        str(bac_multi_index),
        "-l",
        str(bac_paths["lengths"]),
        "-o",
        str(prefix),
    )
    lines = (tmp_out / "extract_x.bed").read_text().strip().splitlines()
    assert len(lines) == 2


def test_get_mum_ranges_flanks_past_index(bac_multi_index):
    idx = sutils.parse_index(str(bac_multi_index), seq_idx=0)
    assert get_mum_ranges_flanks(idx, (9_999_999_999, 9_999_999_999)) is None


def test_extract_no_bounds_fails(
    runner: ShredtoolsRunner,
    bac_sorted,
    bac_multi_index,
    bac_paths,
):
    result = runner.run(
        "extract",
        str(bac_sorted),
        "-s",
        "0",
        "-r",
        "GCF_000009045.1_ASM904v1_genomic.fna:999999999-1000000000",
        "-b",
        str(bac_multi_index),
        "-l",
        str(bac_paths["lengths"]),
        check=False,
    )
    assert result.returncode == 1


def test_extract_plot_requires_output(
    runner: ShredtoolsRunner, bac_sorted, bac_multi_index, bac_paths
):
    result = runner.run(
        "extract",
        str(bac_sorted),
        "-s",
        "0",
        "-r",
        BAC_REGION,
        "-b",
        str(bac_multi_index),
        "-l",
        str(bac_paths["lengths"]),
        "--plot",
        check=False,
    )
    assert result.returncode == 1
    assert "Output prefix required" in result.stderr


@pytest.mark.slow
def test_extract_hg(
    runner: ShredtoolsRunner,
    hg_paths,
    hg_multi_index,
    tmp_out,
):
    prefix = tmp_out / "hg_extract"
    runner.run_ok(
        "extract",
        str(hg_paths["bumbl"]),
        "-s",
        "0",
        "-r",
        "CHM13#0#chr1:1000000-2000000",
        "-b",
        str(hg_multi_index),
        "-l",
        str(hg_paths["lengths"]),
        "-o",
        str(prefix),
    )
    assert (tmp_out / "hg_extract.bed").is_file()
