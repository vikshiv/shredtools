# shredtools tests

Integration and unit tests for all CLI subcommands.

## Quick start

```bash
./scripts/run_tests.sh                  # all tests
./scripts/run_tests.sh -m "not slow"    # bac only (~12s)
./scripts/run_tests.sh -m slow          # hg_test subset (~5s)
./scripts/run_tests.sh -k subset -v     # single module, verbose
```

Requires the `shredtools` conda env and fixtures under `test_data/` (see `test_data/README.md`).

## Layout

| File | Coverage |
|------|----------|
| `test_cli.py` | Help, version, broken pipe |
| `test_stats.py` | `stats` human + JSON |
| `test_sort_filter.py` | `sort`, `filter` |
| `test_index.py` | `index`, checksum, index parsing |
| `test_utils.py` | Coordinate conversion, index helpers, sidecars |
| `test_subset.py` | `subset` happy path + edge cases |
| `test_extract.py` | `extract` BED output + edge cases |
| `test_shred_enhance_fasta.py` | `shred`, `enhance`, `fasta` |

## Markers

- `@pytest.mark.slow` — uses `hg_test.bumbl` (~294 MB)
- `@pytest.mark.requires_refs` — needs FASTA paths from `.lengths` to exist locally

## Environment variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `SHREDTOOLS_TEST_DATA` | `<repo>/test_data` | Fixture directory |
| `SHREDTOOLS_CONDA_ENV` | `shredtools` | Conda env for `run_tests.sh` |
