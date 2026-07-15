"""Tests for shredtools index."""

from __future__ import annotations

import pytest

from shredtools import utils as sutils
from tests.conftest import ShredtoolsRunner


def test_index_unsorted_fails(runner: ShredtoolsRunner, bac_paths, tmp_out):
    out = tmp_out / "bad.bi"
    result = runner.run(
        "index", str(bac_paths["bumbl"]), "--multi", "-o", str(out), check=False
    )
    assert result.returncode == 1
    assert "Sortedness check failed" in result.stderr


def test_index_multi(runner: ShredtoolsRunner, bac_sorted, bac_multi_index):
    assert bac_multi_index.is_file()
    idx = sutils.parse_index(str(bac_multi_index), seq_idx=0)
    assert int(idx.num_seqs) == 5
    assert int(idx.bin_width) == 1_000_000


def test_index_single(runner: ShredtoolsRunner, bac_sorted, bac_single_index):
    idx = sutils.parse_index(str(bac_single_index), seq_idx=0)
    assert int(idx.seq_idx) == 0
    assert idx.offsets.size >= 2


def test_index_single_seq_idx_mismatch(runner: ShredtoolsRunner, bac_single_index):
    with pytest.raises(AssertionError, match="seq_idx mismatch"):
        sutils.parse_index(str(bac_single_index), seq_idx=1)


def test_verify_bumbl_index(runner: ShredtoolsRunner, bac_sorted, bac_multi_index):
    assert sutils.verify_bumbl_index(str(bac_sorted), str(bac_multi_index)) is True


def test_index_custom_bin_width(runner: ShredtoolsRunner, bac_sorted, tmp_out):
    out = tmp_out / "narrow.bi"
    runner.run_ok(
        "index",
        str(bac_sorted),
        "--multi",
        "-w",
        "100000",
        "-o",
        str(out),
    )
    idx = sutils.parse_index(str(out), seq_idx=0)
    assert int(idx.bin_width) == 100_000


def test_checksum_matches_index(bac_sorted, bac_multi_index):
    _n, expected = sutils.checksum_from_bumbl(str(bac_sorted))
    with open(bac_multi_index, "rb") as fin:
        fin.seek(8)
        stored = int.from_bytes(fin.read(8), "little")
    assert stored == int(expected)
