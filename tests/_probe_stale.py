#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
import mumemto.utils as m

bac = REPO / "test_data/bac.bumbl"
sorted_b = Path("/tmp/bac_sorted_stale_test.bumbl")
idx = Path("/tmp/bac_stale_test.bi")
lens = REPO / "test_data/bac.lengths"
region = "GCF_000009045.1_ASM904v1_genomic.fna:100000-200000"

if not sorted_b.exists():
    subprocess.run(["shredtools", "sort", str(bac), "-s", "0", "-o", str(sorted_b)], check=True)
if not idx.exists():
    subprocess.run(["shredtools", "index", str(sorted_b), "--multi", "-o", str(idx)], check=True)

subprocess.run(
    ["shredtools", "subset", str(sorted_b), "-s", "0", "-r", region, "-b", str(idx), "-l", str(lens), "-o", "/tmp/good.bumbl"],
    check=True,
)
subprocess.run(
    ["shredtools", "subset", str(bac), "-s", "0", "-r", region, "-b", str(idx), "-l", str(lens), "-o", "/tmp/bad.bumbl"],
    check=True,
)
g = m.MUMdata("/tmp/good.bumbl", sort=False).num_mums
b = m.MUMdata("/tmp/bad.bumbl", sort=False).num_mums
print(f"good={g} bad={b} stale_mismatch={g != b}")
