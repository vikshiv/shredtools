"""CLI entrypoint tests."""

from __future__ import annotations

import subprocess

import pytest

from tests.conftest import ShredtoolsRunner


def test_top_level_help(runner: ShredtoolsRunner):
    result = runner.run_ok("-h")
    assert "commands:" in result.stdout
    assert "subset" in result.stdout
    assert "extract" in result.stdout


def test_version(runner: ShredtoolsRunner):
    result = runner.run_ok("--version")
    assert "0.1.0" in result.stdout


def test_unknown_command(runner: ShredtoolsRunner):
    result = runner.run("not-a-command", check=False)
    assert result.returncode == 1
    assert "Unknown command" in result.stderr


def test_command_help(runner: ShredtoolsRunner):
    for cmd in ("subset", "extract", "index", "sort", "filter", "stats", "shred", "enhance", "fasta"):
        result = runner.run(cmd, "-h")
        assert result.returncode == 0
        assert "usage:" in result.stdout


def test_broken_pipe_exits_cleanly(runner: ShredtoolsRunner, bac_sorted, bac_multi_index, bac_paths):
    """Piping stdout to head should not print a traceback."""
    cmd = [
        runner.exe,
        "subset",
        str(bac_sorted),
        "-s",
        "0",
        "-r",
        "GCF_000009045.1_ASM904v1_genomic.fna:100000-200000",
        "-b",
        str(bac_multi_index),
        "-l",
        str(bac_paths["lengths"]),
    ]
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert proc.stdout is not None
    proc.stdout.read(64)
    proc.stdout.close()
    proc.wait(timeout=30)
    stderr = proc.stderr.read() if proc.stderr else ""
    # 0 = ok, 141 = SIGPIPE (128+13), -13 = SIGPIPE on some platforms
    assert proc.returncode in (0, 141, -13), stderr
