"""Tests for shred, enhance, and fasta commands."""

from __future__ import annotations

import pytest

from tests.conftest import ShredtoolsRunner, refs_available


def test_shred_creates_beds(
    runner: ShredtoolsRunner,
    bac_sorted,
    bac_paths,
    tmp_out,
):
    out_dir = tmp_out / "shreds"
    runner.run_ok(
        "shred",
        str(bac_sorted),
        "-l",
        str(bac_paths["lengths"]),
        "-o",
        str(out_dir),
        "--shred-size",
        "1000000",
    )
    beds = list(out_dir.glob("*.bed"))
    assert len(beds) >= 2


def test_enhance_no_gaps_fast(
    runner: ShredtoolsRunner,
    bac_sorted,
    bac_paths,
    tmp_out,
):
    """Large min_gap_length => zero gaps enhanced; exercises CLI without mumemto work."""
    out = tmp_out / "enhanced.mums"
    result = runner.run_ok(
        "enhance",
        str(bac_sorted),
        "999999999",
        "-o",
        str(out),
        "-l",
        str(bac_paths["lengths"]),
    )
    assert out.is_file()
    assert "Gaps to enhance: 0" in result.stderr


def test_enhance_help_outlier_and_no_partial(runner: ShredtoolsRunner):
    result = runner.run_ok("enhance", "-h")
    assert "--no-partial" in result.stdout
    assert "--flanks" in result.stdout
    assert "suggested value: 5" in result.stdout
    assert "Default: 5" not in result.stdout


def test_enhance_no_partial_conflicts_with_outlier(tmp_path):
    from shredtools.enhance import parse_arguments

    mum = tmp_path / "x.mums"
    mum.write_text("")
    lens = tmp_path / "x.lengths"
    lens.write_text("")
    with pytest.raises(SystemExit) as exc:
        parse_arguments(
            [
                str(mum),
                "20",
                "-l",
                str(lens),
                "--outlier-factor",
                "5",
                "--no-partial",
            ]
        )
    assert exc.value.code == 1


def test_enhance_gap_no_partial_skips_incomplete():
    from shredtools.enhance import enhance_gap

    # One assembly missing → would be a partial MUM; --no-partial drops it.
    gap_seqs = [(0, "ACGTACGTACGT"), None, (10, "ACGTACGTACGT")]
    assert enhance_gap(gap_seqs, min_match_len=4, no_partial=True) is None


def test_flank_coords_forward_and_reverse():
    """Terminal flanks use orientation-aware first/last block endpoints."""
    import numpy as np
    import mumemto.utils as mutils
    from shredtools.enhance import compute_gap_coords_and_lengths

    # One block of two MUMs across 3 genomes; genome 1 is reverse-oriented.
    lengths = np.array([10, 10], dtype=np.uint32)
    starts = np.array(
        [
            [100, 200, 50],
            [200, 50, 150],
        ],
        dtype=np.int64,
    )
    strands = np.array(
        [
            [True, False, True],
            [True, False, True],
        ],
        dtype=bool,
    )
    mums = mutils.MUMdata.from_arrays(lengths, starts, strands, blocks=[(0, 1)])
    total_lens = [1000, 1000, 1000]

    gaps_plain, coords_plain, lens_plain = compute_gap_coords_and_lengths(mums)
    assert gaps_plain == [(0, 1)]
    assert coords_plain.shape == (1, 3, 2)

    gaps, coords, gap_lengths = compute_gap_coords_and_lengths(
        mums, include_flanks=True, total_lens=total_lens
    )
    assert len(gaps) == 3  # 1 internal + left + right (no inter with 1 block)
    assert gaps[0] == (0, 1)
    assert gaps[1] == ("left", 0, 1)
    assert gaps[2] == ("right", 0, 1)
    np.testing.assert_array_equal(coords[0], coords_plain[0])
    np.testing.assert_array_equal(gap_lengths[0], lens_plain[0])

    # Forward genomes: left [0, first_start), right [last_end, genome_len)
    np.testing.assert_array_equal(coords[1, 0], [0, 100])
    np.testing.assert_array_equal(coords[2, 0], [210, 1000])
    np.testing.assert_array_equal(coords[1, 2], [0, 50])
    np.testing.assert_array_equal(coords[2, 2], [160, 1000])

    # Reverse genome: leftmost is last MUM start; rightmost is first MUM end
    np.testing.assert_array_equal(coords[1, 1], [0, 50])
    np.testing.assert_array_equal(coords[2, 1], [210, 1000])
    np.testing.assert_array_equal(gap_lengths[1], [100, 50, 50])
    np.testing.assert_array_equal(gap_lengths[2], [790, 790, 840])


def test_flank_coords_singleton_block():
    """Singleton blocks have no internal gaps but still emit terminal flanks."""
    import numpy as np
    import mumemto.utils as mutils
    from shredtools.enhance import compute_gap_coords_and_lengths

    lengths = np.array([20], dtype=np.uint32)
    starts = np.array([[10, 30]], dtype=np.int64)
    strands = np.ones((1, 2), dtype=bool)
    mums = mutils.MUMdata.from_arrays(lengths, starts, strands, blocks=[(0, 0)])
    total_lens = [100, 200]

    with pytest.raises(SystemExit):
        compute_gap_coords_and_lengths(mums)

    gaps, coords, gap_lengths = compute_gap_coords_and_lengths(
        mums, include_flanks=True, total_lens=total_lens
    )
    assert gaps == [("left", 0, 0), ("right", 0, 0)]
    np.testing.assert_array_equal(coords[0], [[0, 10], [0, 30]])
    np.testing.assert_array_equal(coords[1], [[30, 100], [50, 200]])
    np.testing.assert_array_equal(gap_lengths[0], [10, 30])
    np.testing.assert_array_equal(gap_lengths[1], [70, 150])


def test_flank_coords_multi_block_terminal_and_inter():
    """Multi-block: only outer terminals + inter-block gaps (not per-block to 0/end)."""
    import numpy as np
    import mumemto.utils as mutils
    from shredtools.enhance import compute_gap_coords_and_lengths

    # 3 blocks: (0,1), (2,2) singleton, (3,4). Two genomes, all forward.
    lengths = np.array([10, 10, 5, 10, 10], dtype=np.uint32)
    starts = np.array(
        [
            [100, 100],
            [200, 200],
            [400, 400],  # singleton block
            [600, 600],
            [700, 700],
        ],
        dtype=np.int64,
    )
    strands = np.ones((5, 2), dtype=bool)
    mums = mutils.MUMdata.from_arrays(
        lengths, starts, strands, blocks=[(0, 1), (2, 2), (3, 4)]
    )
    total_lens = [1000, 1000]

    gaps_plain, coords_plain, _ = compute_gap_coords_and_lengths(mums)
    assert gaps_plain == [(0, 1), (3, 4)]  # internals only; singleton has none

    gaps, coords, gap_lengths = compute_gap_coords_and_lengths(
        mums, include_flanks=True, total_lens=total_lens
    )
    # 2 internals + left + 2 inter + right
    assert gaps == [
        (0, 1),
        (3, 4),
        ("left", 0, 1),
        ("inter", 1, 2),
        ("inter", 2, 3),
        ("right", 3, 4),
    ]
    np.testing.assert_array_equal(coords[:2], coords_plain)

    # Left terminal from first block only: [0, 100)
    np.testing.assert_array_equal(coords[2], [[0, 100], [0, 100]])
    # Inter block0→block1: after MUM1 end (210) to MUM2 start (400)
    np.testing.assert_array_equal(coords[3], [[210, 400], [210, 400]])
    # Inter block1→block2: after MUM2 end (405) to MUM3 start (600)
    np.testing.assert_array_equal(coords[4], [[405, 600], [405, 600]])
    # Right terminal from last block only: [710, 1000)
    np.testing.assert_array_equal(coords[5], [[710, 1000], [710, 1000]])

    # Intermediate blocks must not emit [0, ...) or [..., genome_end) flanks
    kinds = [g[0] if isinstance(g[0], str) else "gap" for g in gaps]
    assert kinds.count("left") == 1
    assert kinds.count("right") == 1
    assert gap_lengths[2, 0] == 100
    assert gap_lengths[5, 0] == 290


def test_flanks_multi_block_no_duplicate_input_mums(tmp_path):
    """Regression: --flanks on multi-block input must not re-emit exact input rows."""
    import numpy as np
    import mumemto.utils as mutils
    from mumemto import mum as run_mumemto
    from shredtools.enhance import main as enhance_main

    # Two collinear blocks with a reciprocal rearrangement so find_coll_blocks
    # yields >=2 blocks (not one merged run). Divergent spacers keep MUMs split.
    s0, s1 = "ACGTACGTACGTACGTAC", "GGAAGGAAGGAAAGGAA"
    s2, s3 = "TGCATGCATGCATGCATG", "CCTTCCTTCCTTCCTTC"
    blk0_a = s0 + "AAAA" + s1
    blk0_b = s0 + "TTTT" + s1
    blk1_a = s2 + "CCCC" + s3
    blk1_b = s2 + "GGGG" + s3
    pad = "GGGCCCCG"
    mid_a = "AAAATTTTAAAATTTT"
    mid_b = "TTTTTAAATTTTTTAA"
    seq_a = pad + blk0_a + mid_a + blk1_a + pad
    seq_b = pad + blk1_b + mid_b + blk0_b + pad

    fa_a = tmp_path / "a.fa"
    fa_b = tmp_path / "b.fa"
    fa_a.write_text(f">a\n{seq_a}\n")
    fa_b.write_text(f">b\n{seq_b}\n")
    n = len(seq_a)
    lengths_path = tmp_path / "toy.lengths"
    lengths_path.write_text(
        f"{fa_a} * {n}\n{fa_a} a {n}\n{fa_b} * {n}\n{fa_b} b {n}\n"
    )

    raw = run_mumemto([[seq_a], [seq_b]], min_match_len=12)
    lengths = np.array([r[0] for r in raw], dtype=np.uint32)
    starts = np.array([r[1] for r in raw], dtype=np.int64)
    strands = np.array([r[2] for r in raw], dtype=bool)
    mums = mutils.MUMdata.from_arrays(lengths, starts, strands)
    mums.sort(ref_col=0, copy=False)
    mums.blocks = mutils.find_coll_blocks(mums)
    assert len(mums.blocks) >= 2, f"need multi-block toy, got {mums.blocks}"

    mum_path = tmp_path / "toy.mums"
    mums.write_mums(str(mum_path))
    out_path = tmp_path / "toy.enhanced.mums"

    enhance_main(
        [
            str(mum_path),
            "0",
            "--flanks",
            "-l",
            str(lengths_path),
            "-o",
            str(out_path),
            "-t",
            "1",
        ]
    )

    out = mutils.MUMdata(str(out_path), sort=False)
    inp = mutils.MUMdata(str(mum_path), sort=False)

    def keys(m):
        return [
            (int(m.lengths[i]), tuple(int(x) for x in m.starts[i]))
            for i in range(len(m))
        ]

    in_keys = keys(inp)
    out_keys = keys(out)
    assert len(in_keys) == len(set(in_keys))
    for k in in_keys:
        assert out_keys.count(k) == 1, f"input MUM duplicated in output: {k}"


def test_parse_arguments_flanks_default_off(tmp_path):
    from shredtools.enhance import parse_arguments

    mum = tmp_path / "x.mums"
    mum.write_text("")
    lens = tmp_path / "x.lengths"
    lens.write_text("")
    args = parse_arguments([str(mum), "20", "-l", str(lens)])
    assert args.flanks is False
    args = parse_arguments([str(mum), "20", "-l", str(lens), "--flanks"])
    assert args.flanks is True


@pytest.mark.requires_refs
def test_fasta_from_extract_bed(
    runner: ShredtoolsRunner,
    bac_sorted,
    bac_multi_index,
    bac_paths,
    tmp_out,
):
    if not refs_available(bac_paths["lengths"]):
        pytest.skip("Reference FASTAs from bac.lengths not available on this machine")

    prefix = tmp_out / "regions"
    runner.run_ok(
        "extract",
        str(bac_sorted),
        "-s",
        "0",
        "-r",
        "GCF_000009045.1_ASM904v1_genomic.fna:100000-200000",
        "-b",
        str(bac_multi_index),
        "-l",
        str(bac_paths["lengths"]),
        "-o",
        str(prefix),
    )
    fasta_dir = tmp_out / "fasta_out"
    runner.run_ok("fasta", str(tmp_out / "regions.bed"), "-o", str(fasta_dir))
    fastas = list(fasta_dir.glob("*.fa"))
    assert len(fastas) >= 1
    assert fastas[0].stat().st_size > 0


def test_fasta_missing_bed(runner: ShredtoolsRunner, tmp_out):
    result = runner.run("fasta", str(tmp_out / "missing.bed"), "-o", str(tmp_out / "out"), check=False)
    assert result.returncode == 1
