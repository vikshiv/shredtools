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


def test_filter_help_mentions_partial(runner: ShredtoolsRunner):
    result = runner.run_ok("filter", "-h")
    assert "--partial" in result.stdout
    assert "collinear" in result.stdout.lower()
    assert "partial" in result.stdout.lower()


def test_filter_partial_drops_absent_starts(tmp_path):
    import numpy as np
    import mumemto.utils as mutils
    from shredtools.filter_collinear import filter_partial_mums

    lengths = np.array([10, 12, 8], dtype=np.uint32)
    starts = np.array(
        [
            [0, 0, 0],
            [20, -1, 20],  # partial
            [40, 40, 40],
        ],
        dtype=np.int64,
    )
    strands = np.ones((3, 3), dtype=bool)
    mums = mutils.MUMdata.from_arrays(lengths, starts, strands)
    out = filter_partial_mums(mums)
    assert len(out) == 2
    assert not np.any(out.starts == -1)


def test_filter_partial_inplace_preserves_mums_text(tmp_path):
    import numpy as np
    import mumemto.utils as mutils
    from shredtools.filter_collinear import main

    lengths = np.array([10, 12, 8], dtype=np.uint32)
    starts = np.array(
        [
            [0, 0, 0],
            [20, -1, 20],
            [40, 40, 40],
        ],
        dtype=np.int64,
    )
    strands = np.ones((3, 3), dtype=bool)
    path = tmp_path / "plus.enhanced.mums"
    mutils.MUMdata.from_arrays(lengths, starts, strands).write_mums(str(path))

    main([str(path), "--partial", "-i"])

    text = path.read_text()
    assert "\t" in text  # text .mums, not binary bumbl
    assert "-1" not in text
    reloaded = mutils.MUMdata(str(path), sort=False)
    assert len(reloaded) == 2
