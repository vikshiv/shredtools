"""Unit tests for shredtools.utils helpers."""

from __future__ import annotations

import pytest

import mumemto.utils as mutils
from shredtools import utils as sutils


BAC_REGION = "GCF_000009045.1_ASM904v1_genomic.fna:100000-200000"
BAC_CONTIG = "GCF_000009045.1_ASM904v1_genomic.fna"


def test_is_url():
    assert sutils.is_url("https://example.com/a.bumbl")
    assert sutils.is_url("http://x/y")
    assert not sutils.is_url("/local/path.bumbl")
    assert not sutils.is_url("ftp://x")


def test_convert_local_to_global_coords(bac_paths):
    names = mutils.get_contig_names(str(bac_paths["lengths"]))
    lengths = mutils.get_sequence_lengths(str(bac_paths["lengths"]), multilengths=True)[0]
    start, end = sutils.convert_local_to_global_coords(
        BAC_REGION, names[0], lengths
    )
    assert start == 100_000
    assert end == 199_999  # inclusive global end; region is half-open [start, end)


def test_convert_local_invalid_contig(bac_paths):
    names = mutils.get_contig_names(str(bac_paths["lengths"]))
    lengths = mutils.get_sequence_lengths(str(bac_paths["lengths"]), multilengths=True)[0]
    with pytest.raises(AssertionError):
        sutils.convert_local_to_global_coords("not_a_contig:1-100", names[0], lengths)


def test_convert_local_out_of_range(bac_paths):
    names = mutils.get_contig_names(str(bac_paths["lengths"]))
    lengths = mutils.get_sequence_lengths(str(bac_paths["lengths"]), multilengths=True)[0]
    with pytest.raises(AssertionError):
        sutils.convert_local_to_global_coords(
            f"{BAC_CONTIG}:999999999-1000000000", names[0], lengths
        )


def test_get_mum_ranges_region(bac_multi_index, bac_paths):
    names = mutils.get_contig_names(str(bac_paths["lengths"]))
    lengths = mutils.get_sequence_lengths(str(bac_paths["lengths"]), multilengths=True)[0]
    coords = sutils.convert_local_to_global_coords(BAC_REGION, names[0], lengths)
    idx = sutils.parse_index(str(bac_multi_index), seq_idx=0)
    ranges = sutils.get_mum_ranges_region(idx, coords)
    assert ranges is not None
    assert len(ranges) >= 1


def test_get_mum_ranges_region_past_index(bac_multi_index):
    idx = sutils.parse_index(str(bac_multi_index), seq_idx=0)
    # Coordinate far beyond indexed bins
    ranges = sutils.get_mum_ranges_region(idx, (9_999_999_999, 9_999_999_999))
    assert ranges is None


def test_parse_bumbl_range_roundtrip(bac_sorted, bac_multi_index, bac_paths):
    names = mutils.get_contig_names(str(bac_paths["lengths"]))
    lengths = mutils.get_sequence_lengths(str(bac_paths["lengths"]), multilengths=True)[0]
    coords = sutils.convert_local_to_global_coords(BAC_REGION, names[0], lengths)
    idx = sutils.parse_index(str(bac_multi_index), seq_idx=0)
    ranges = sutils.get_mum_ranges_region(idx, coords)
    partial = sutils.parse_bumbl_range(str(bac_sorted), ranges)
    full = mutils.MUMdata(str(bac_sorted), sort=False)
    assert partial.num_seqs == full.num_seqs
    assert 0 < partial.num_mums <= full.num_mums


def test_resolve_bumbl_sidecars_with_explicit_paths(bac_paths, bac_multi_index):
    bi, lens = sutils.resolve_bumbl_sidecars(
        str(bac_paths["bumbl"]),
        bi=str(bac_multi_index),
        lens=str(bac_paths["lengths"]),
    )
    assert bi == str(bac_multi_index)
    assert lens == str(bac_paths["lengths"])


def test_resolve_bumbl_sidecars_missing(tmp_path):
    missing = tmp_path / "nope.bumbl"
    with pytest.raises(SystemExit):
        sutils.resolve_bumbl_sidecars(str(missing))
