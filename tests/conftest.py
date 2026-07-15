"""Shared fixtures and helpers for shredtools integration tests."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TEST_DATA = REPO_ROOT / "test_data"


def test_data_dir() -> Path:
    return Path(os.environ.get("SHREDTOOLS_TEST_DATA", DEFAULT_TEST_DATA))


def require_file(path: Path, label: str) -> Path:
    if not path.is_file():
        pytest.skip(f"{label} not found: {path}")
    return path


class ShredtoolsRunner:
    """Invoke the shredtools CLI as a subprocess."""

    def __init__(self, exe: str = "shredtools"):
        self.exe = exe

    def run(
        self,
        *args: str,
        check: bool = True,
        text: bool = True,
        capture: bool = True,
        input: str | None = None,
    ) -> subprocess.CompletedProcess[str]:
        cmd = [self.exe, *args]
        return subprocess.run(
            cmd,
            check=check,
            text=text,
            capture_output=capture,
            input=input,
        )

    def run_ok(self, *args: str, **kwargs) -> subprocess.CompletedProcess[str]:
        result = self.run(*args, **kwargs)
        assert result.returncode == 0, (
            f"command failed ({result.returncode}): shredtools {' '.join(args)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
        return result


@pytest.fixture(scope="session")
def runner() -> ShredtoolsRunner:
    return ShredtoolsRunner()


@pytest.fixture(scope="session")
def data_dir() -> Path:
    d = test_data_dir()
    if not d.is_dir():
        pytest.skip(f"Test data directory not found: {d}")
    return d


@pytest.fixture(scope="session")
def bac_paths(data_dir: Path) -> dict[str, Path]:
    return {
        "bumbl": require_file(data_dir / "bac.bumbl", "bac.bumbl"),
        "lengths": require_file(data_dir / "bac.lengths", "bac.lengths"),
    }


@pytest.fixture(scope="session")
def hg_paths(data_dir: Path) -> dict[str, Path]:
    return {
        "bumbl": require_file(data_dir / "hg_test.bumbl", "hg_test.bumbl"),
        "lengths": require_file(data_dir / "hg_test.lengths", "hg_test.lengths"),
    }


@pytest.fixture(scope="session")
def work_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    d = tmp_path_factory.mktemp("shredtools_session")
    return Path(d)


@pytest.fixture(scope="session")
def bac_sorted(work_dir: Path, bac_paths: dict[str, Path], runner: ShredtoolsRunner) -> Path:
    """Sorted copy of bac.bumbl (index requires sorted rows)."""
    out = work_dir / "bac_sorted.bumbl"
    if not out.exists():
        runner.run_ok("sort", str(bac_paths["bumbl"]), "-s", "0", "-o", str(out))
    return out


@pytest.fixture(scope="session")
def bac_multi_index(work_dir: Path, bac_sorted: Path, runner: ShredtoolsRunner) -> Path:
    out = work_dir / "bac_multi.bi"
    if not out.exists():
        runner.run_ok("index", str(bac_sorted), "--multi", "-o", str(out))
    return out


@pytest.fixture(scope="session")
def bac_single_index(work_dir: Path, bac_sorted: Path, runner: ShredtoolsRunner) -> Path:
    out = work_dir / "bac_single.bi"
    if not out.exists():
        runner.run_ok(
            "index", str(bac_sorted), "--single", "-s", "0", "-o", str(out)
        )
    return out


@pytest.fixture(scope="session")
def hg_multi_index(work_dir: Path, hg_paths: dict[str, Path], runner: ShredtoolsRunner) -> Path:
    out = work_dir / "hg_multi.bi"
    if not out.exists():
        runner.run_ok("index", str(hg_paths["bumbl"]), "--multi", "-o", str(out))
    return out


@pytest.fixture
def tmp_out(tmp_path: Path) -> Path:
    return tmp_path


def load_mum_header(path: Path) -> dict:
    from shredtools.utils import read_mum_header

    info = read_mum_header(str(path))
    return {
        "n_mums": info.n_mums,
        "n_seqs": info.n_seqs,
        "blocks_stored": info.blocks_stored,
    }


def refs_available(lengths_path: Path) -> bool:
    """True if every FASTA path in a multilengths file exists."""
    import mumemto.utils as mutils

    try:
        paths = mutils.get_seq_paths(str(lengths_path))
    except Exception:
        return False
    return all(Path(p).is_file() for p in paths)


@pytest.fixture(scope="session")
def bac_refs_ok(bac_paths: dict[str, Path]) -> bool:
    return refs_available(bac_paths["lengths"])


def parse_stats_json(runner: ShredtoolsRunner, bumbl: Path) -> dict:
    result = runner.run_ok("stats", str(bumbl), "--json")
    return json.loads(result.stdout)
