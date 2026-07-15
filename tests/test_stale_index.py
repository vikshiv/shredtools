"""Stale .bumbl.bi index produces wrong subset rows (no runtime checksum check)."""

from __future__ import annotations

import subprocess

import pytest
import mumemto.utils as mutils

from tests.conftest import get_test_data_dir

BAC_REGION = "GCF_000009045.1_ASM904v1_genomic.fna:100000-200000"


@pytest.fixture(scope="module")
def stale_index_setup(tmp_path_factory):
    """Index built from sorted bac; query uses unsorted original."""
    data = get_test_data_dir()
    bac = data / "bac.bumbl"
    lens = data / "bac.lengths"
    if not bac.is_file():
        pytest.skip("bac.bumbl not available")

    work = tmp_path_factory.mktemp("stale_index")
    sorted_b = work / "bac_sorted.bumbl"
    index = work / "bac.bi"
    good_out = work / "good.bumbl"
    bad_out = work / "bad.bumbl"

    subprocess.run(
        ["shredtools", "sort", str(bac), "-s", "0", "-o", str(sorted_b)],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        ["shredtools", "index", str(sorted_b), "--multi", "-o", str(index)],
        check=True,
        capture_output=True,
    )
    region_args = [
        "-s",
        "0",
        "-r",
        BAC_REGION,
        "-b",
        str(index),
        "-l",
        str(lens),
    ]
    subprocess.run(
        ["shredtools", "subset", str(sorted_b), *region_args, "-o", str(good_out)],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        ["shredtools", "subset", str(bac), *region_args, "-o", str(bad_out)],
        check=True,
        capture_output=True,
    )
    return {
        "sorted_b": sorted_b,
        "index": index,
        "lens": lens,
        "bac": bac,
        "good_out": good_out,
        "bad_out": bad_out,
    }


def test_stale_index_subset_row_count_mismatch(stale_index_setup):
    """Document current behavior: wrong index silently returns incorrect row counts."""
    good = mutils.MUMdata(str(stale_index_setup["good_out"]), sort=False).num_mums
    bad = mutils.MUMdata(str(stale_index_setup["bad_out"]), sort=False).num_mums
    assert good == 519  # calibrated on bac test data
    assert bad != good
    assert bad < good


def test_stale_index_ground_truth_matches_sorted_file(stale_index_setup):
    """Sorted bumbl + matching index should match in-memory overlap filter."""
    from shredtools import utils as sutils

    setup = stale_index_setup
    subset = mutils.MUMdata(str(setup["good_out"]), sort=False)
    names = mutils.get_contig_names(str(setup["lens"]))
    lengths = mutils.get_sequence_lengths(str(setup["lens"]), multilengths=True)[0]
    coords = sutils.convert_local_to_global_coords(BAC_REGION, names[0], lengths)
    region_start, region_end_incl = coords
    region_end_excl = region_end_incl + 1
    full = mutils.MUMdata(str(setup["sorted_b"]), sort=False)
    starts = full.starts[:, 0]
    overlap = (starts >= 0) & (starts < region_end_excl) & (
        starts + full.lengths > region_start
    )
    assert subset.num_mums == int(overlap.sum())
