"""Tests for shredtools subset."""

from __future__ import annotations

import pytest

import mumemto.utils as mutils
from tests.conftest import ShredtoolsRunner, load_mum_header

BAC_REGION = "GCF_000009045.1_ASM904v1_genomic.fna:100000-200000"


def test_subset_to_bumbl(
    runner: ShredtoolsRunner,
    bac_sorted,
    bac_multi_index,
    bac_paths,
    tmp_out,
):
    out = tmp_out / "subset.bumbl"
    runner.run_ok(
        "subset",
        str(bac_sorted),
        "-s",
        "0",
        "-r",
        BAC_REGION,
        "-b",
        str(bac_multi_index),
        "-l",
        str(bac_paths["lengths"]),
        "-v",
        "-o",
        str(out),
    )
    header = load_mum_header(out)
    assert 0 < header["n_mums"] < 6848
    assert header["n_seqs"] == 5


def test_subset_stdout_mums(
    runner: ShredtoolsRunner,
    bac_sorted,
    bac_multi_index,
    bac_paths,
):
    result = runner.run_ok(
        "subset",
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
    assert len(result.stdout) > 50
    assert "kept" not in result.stderr


def test_subset_verbose_stderr(
    runner: ShredtoolsRunner,
    bac_sorted,
    bac_multi_index,
    bac_paths,
    tmp_out,
):
    out = tmp_out / "subset.bumbl"
    result = runner.run_ok(
        "subset",
        str(bac_sorted),
        "-s",
        "0",
        "-r",
        BAC_REGION,
        "-b",
        str(bac_multi_index),
        "-l",
        str(bac_paths["lengths"]),
        "-v",
        "-o",
        str(out),
    )
    assert "loaded" in result.stderr
    assert "kept" in result.stderr


def test_subset_matches_manual_filter(
    runner: ShredtoolsRunner,
    bac_sorted,
    bac_multi_index,
    bac_paths,
    tmp_out,
):
    """Subset row count should match an in-memory overlap filter."""
    from shredtools import utils as sutils

    out = tmp_out / "subset.bumbl"
    runner.run_ok(
        "subset",
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
        str(out),
    )
    subset = mutils.MUMdata(str(out), sort=False)
    names = mutils.get_contig_names(str(bac_paths["lengths"]))
    lengths = mutils.get_sequence_lengths(str(bac_paths["lengths"]), multilengths=True)[0]
    coords = sutils.convert_local_to_global_coords(BAC_REGION, names[0], lengths)
    region_start, region_end_incl = coords
    region_end_excl = region_end_incl + 1
    full = mutils.MUMdata(str(bac_sorted), sort=False)
    starts = full.starts[:, 0]
    overlap = (starts >= 0) & (starts < region_end_excl) & (
        starts + full.lengths > region_start
    )
    assert subset.num_mums == int(overlap.sum())


def test_subset_invalid_seq_idx(
    runner: ShredtoolsRunner,
    bac_sorted,
    bac_multi_index,
    bac_paths,
):
    result = runner.run(
        "subset",
        str(bac_sorted),
        "-s",
        "99",
        "-r",
        BAC_REGION,
        "-b",
        str(bac_multi_index),
        "-l",
        str(bac_paths["lengths"]),
        check=False,
    )
    assert result.returncode == 1
    assert "invalid" in result.stderr.lower()


def test_subset_invalid_region(
    runner: ShredtoolsRunner,
    bac_sorted,
    bac_multi_index,
    bac_paths,
):
    result = runner.run(
        "subset",
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


def test_subset_single_index(
    runner: ShredtoolsRunner,
    bac_sorted,
    bac_single_index,
    bac_paths,
    tmp_out,
):
    out = tmp_out / "subset_single.bumbl"
    runner.run_ok(
        "subset",
        str(bac_sorted),
        "-s",
        "0",
        "-r",
        BAC_REGION,
        "-b",
        str(bac_single_index),
        "-l",
        str(bac_paths["lengths"]),
        "-o",
        str(out),
    )
    assert load_mum_header(out)["n_mums"] > 0


@pytest.mark.slow
def test_subset_hg(
    runner: ShredtoolsRunner,
    hg_paths,
    hg_multi_index,
    tmp_out,
):
    out = tmp_out / "hg_subset.bumbl"
    runner.run_ok(
        "subset",
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
        str(out),
    )
    header = load_mum_header(out)
    assert header["n_mums"] > 0
    assert header["n_seqs"] == 6
