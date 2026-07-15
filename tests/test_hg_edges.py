"""Human-genome edge-case tests derived from hg_test.lengths coordinate analysis."""

from __future__ import annotations

import pytest

import mumemto.utils as mutils
from shredtools import utils as sutils
from shredtools.extract_from_mums import get_mum_ranges_flanks
from tests.conftest import ShredtoolsRunner, load_mum_header
from tests.hg_edge_regions import (
    CHR1_MID,
    CHR1_SPARSE_TAIL,
    CHR1_START,
    CHRM_EMPTY,
    CHRY_TAIL,
)

pytestmark = pytest.mark.slow


def _coords(region: str, hg_paths) -> tuple[int, int]:
    names = mutils.get_contig_names(str(hg_paths["lengths"]))
    lengths = mutils.get_sequence_lengths(str(hg_paths["lengths"]), multilengths=True)[0]
    return sutils.convert_local_to_global_coords(region, names[0], lengths)


def test_hg_chrM_past_index_bins(hg_paths, hg_multi_index):
    """chrM sits at the genome tail; global coords exceed max_bin."""
    coords = _coords(CHRM_EMPTY, hg_paths)
    idx = sutils.parse_index(str(hg_multi_index), seq_idx=0)
    assert coords[0] > 3_117_000_000
    assert sutils.get_mum_ranges_region(idx, coords) is None
    assert get_mum_ranges_flanks(idx, coords) is None


def test_hg_subset_empty_chrM_exits_zero(
    runner: ShredtoolsRunner, hg_paths, hg_multi_index, tmp_out
):
    out = tmp_out / "chrm.bumbl"
    result = runner.run_ok(
        "subset",
        str(hg_paths["bumbl"]),
        "-s",
        "0",
        "-r",
        CHRM_EMPTY,
        "-b",
        str(hg_multi_index),
        "-l",
        str(hg_paths["lengths"]),
        "-v",
        "-o",
        str(out),
    )
    assert "No indexed bins found" in result.stderr
    assert "kept 0 MUM rows" in result.stderr
    assert load_mum_header(out)["n_mums"] == 0


def test_hg_extract_fails_chrM(
    runner: ShredtoolsRunner, hg_paths, hg_multi_index
):
    result = runner.run(
        "extract",
        str(hg_paths["bumbl"]),
        "-s",
        "0",
        "-r",
        CHRM_EMPTY,
        "-b",
        str(hg_multi_index),
        "-l",
        str(hg_paths["lengths"]),
        check=False,
    )
    assert result.returncode == 1
    assert "No bounding MUMs found" in result.stderr


def test_hg_subset_vs_extract_asymmetric_genome_tail(
    runner: ShredtoolsRunner, hg_paths, hg_multi_index, tmp_out
):
    """Same region: subset exit 0 with empty output; extract exit 1."""
    out = tmp_out / "chry.bumbl"
    sub = runner.run_ok(
        "subset",
        str(hg_paths["bumbl"]),
        "-s",
        "0",
        "-r",
        CHRY_TAIL,
        "-b",
        str(hg_multi_index),
        "-l",
        str(hg_paths["lengths"]),
        "-o",
        str(out),
    )
    ext = runner.run(
        "extract",
        str(hg_paths["bumbl"]),
        "-s",
        "0",
        "-r",
        CHRY_TAIL,
        "-b",
        str(hg_multi_index),
        "-l",
        str(hg_paths["lengths"]),
        check=False,
    )
    assert "No indexed bins found" in sub.stderr
    assert load_mum_header(out)["n_mums"] == 0
    assert ext.returncode == 1
    assert "No bounding MUMs found" in ext.stderr


def test_hg_chr1_start_index_bins_but_flanks_fail(hg_paths, hg_multi_index):
    """Telomere gap: bins exist but get_mum_ranges_flanks returns None."""
    coords = _coords(CHR1_START, hg_paths)
    idx = sutils.parse_index(str(hg_multi_index), seq_idx=0)
    assert sutils.get_mum_ranges_region(idx, coords) is not None
    assert get_mum_ranges_flanks(idx, coords) is None


def test_hg_subset_keeps_rows_extract_fails_chr1_start(
    runner: ShredtoolsRunner, hg_paths, hg_multi_index, tmp_out
):
    """chr1:1-50000 — subset finds 6 overlapping rows; extract cannot bound region."""
    out = tmp_out / "chr1_start.bumbl"
    sub = runner.run_ok(
        "subset",
        str(hg_paths["bumbl"]),
        "-s",
        "0",
        "-r",
        CHR1_START,
        "-b",
        str(hg_multi_index),
        "-l",
        str(hg_paths["lengths"]),
        "-v",
        "-o",
        str(out),
    )
    ext = runner.run(
        "extract",
        str(hg_paths["bumbl"]),
        "-s",
        "0",
        "-r",
        CHR1_START,
        "-b",
        str(hg_multi_index),
        "-l",
        str(hg_paths["lengths"]),
        check=False,
    )
    n = load_mum_header(out)["n_mums"]
    assert n > 0, "expected small number of telomere MUMs"
    assert "loaded" in sub.stderr
    assert ext.returncode == 1
    assert "No bounding MUMs found" in ext.stderr


def test_hg_chr1_sparse_tail_subset_loads_but_filters_zero(
    runner: ShredtoolsRunner, hg_paths, hg_multi_index, tmp_out
):
    """Past last chr1 MUM locally: index returns candidates, overlap filter keeps 0."""
    out = tmp_out / "chr1_tail.bumbl"
    result = runner.run_ok(
        "subset",
        str(hg_paths["bumbl"]),
        "-s",
        "0",
        "-r",
        CHR1_SPARSE_TAIL,
        "-b",
        str(hg_multi_index),
        "-l",
        str(hg_paths["lengths"]),
        "-v",
        "-o",
        str(out),
    )
    assert "loaded" in result.stderr
    assert "kept 0 MUM rows" in result.stderr
    assert load_mum_header(out)["n_mums"] == 0


def test_hg_chr1_sparse_tail_extract_empty_bed(
    runner: ShredtoolsRunner, hg_paths, hg_multi_index, tmp_out
):
    """Extract finds margins but all BED rows skipped (cross-contig on every haplotype)."""
    prefix = tmp_out / "chr1_tail"
    result = runner.run_ok(
        "extract",
        str(hg_paths["bumbl"]),
        "-s",
        "0",
        "-r",
        CHR1_SPARSE_TAIL,
        "-b",
        str(hg_multi_index),
        "-l",
        str(hg_paths["lengths"]),
        "-o",
        str(prefix),
    )
    bed = tmp_out / "chr1_tail.bed"
    assert bed.is_file()
    assert bed.read_text().strip() == ""
    assert "left margin" in result.stderr
    assert "Skipping BED line" in result.stderr


def test_hg_chr1_mid_sanity(
    runner: ShredtoolsRunner, hg_paths, hg_multi_index, tmp_out
):
    """Interior region should succeed for both subset and extract."""
    sub_out = tmp_out / "mid.bumbl"
    runner.run_ok(
        "subset",
        str(hg_paths["bumbl"]),
        "-s",
        "0",
        "-r",
        CHR1_MID,
        "-b",
        str(hg_multi_index),
        "-l",
        str(hg_paths["lengths"]),
        "-o",
        str(sub_out),
    )
    assert load_mum_header(sub_out)["n_mums"] > 100

    prefix = tmp_out / "mid_extract"
    runner.run_ok(
        "extract",
        str(hg_paths["bumbl"]),
        "-s",
        "0",
        "-r",
        CHR1_MID,
        "-b",
        str(hg_multi_index),
        "-l",
        str(hg_paths["lengths"]),
        "-o",
        str(prefix),
    )
    bed_lines = (tmp_out / "mid_extract.bed").read_text().strip().splitlines()
    assert len(bed_lines) == 6
