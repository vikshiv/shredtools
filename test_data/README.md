# Test data for shredtools

Integration tests read fixtures from this directory. **Do not commit large binary files** to git.

## Expected files

| File | Size (approx.) | Role |
|------|----------------|------|
| `bac.bumbl` | ~300 KB | Small bacterial pangenome (5 seqs, ~7k MUMs). Primary fast-test dataset. |
| `bac.lengths` | ~1 KB | mumemto multilengths manifest for `bac.bumbl` |
| `hg_test.bumbl` | ~294 MB | Human pangenome subset (6 seqs, ~5.8M MUMs). Marked `@pytest.mark.slow`. |
| `hg_test.lengths` | ~52 KB | Manifest for `hg_test.bumbl` |
| `hg_test.athresh` | ~5.9 GB | Optional sidecar probed by `stats`; not required for most tests |

## Setup

Place your `.bumbl` and `.lengths` files here (or set `SHREDTOOLS_TEST_DATA` to another directory).

Reference FASTAs listed in the `.lengths` files must exist for `fasta` / `enhance` tests that touch real sequence. Those tests skip automatically when paths are missing.

## Running tests

```bash
# From repo root, using the shredtools conda env:
./scripts/run_tests.sh

# Fast tests only (bac dataset):
./scripts/run_tests.sh -m "not slow"

# Include human subset:
./scripts/run_tests.sh -m slow
```

## Commit guidance

| Commit to git | Keep private / gitignore |
|---------------|--------------------------|
| All files under `tests/` | `hg_test.bumbl`, `hg_test.athresh` |
| `test_data/README.md`, `test_data/*.lengths` (optional) | Generated `.bi` indexes |
| `pytest.ini`, `scripts/run_tests.sh` | `tests/output/` |
| Optionally `bac.bumbl` (~300 KB) if you want CI without extra setup | Any local FASTA mirrors |

The **test code** should always be committed. **Large binaries** should stay local or live in object storage with download instructions.
