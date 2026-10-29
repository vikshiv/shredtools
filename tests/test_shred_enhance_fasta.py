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
