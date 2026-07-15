#!/usr/bin/env python3
"""Probe hg_test edge regions — run: conda run -n shredtools python tests/_probe_hg_edges.py"""
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import mumemto.utils as m
from shredtools import utils as sutils
from shredtools.extract_from_mums import get_mum_ranges_flanks

lens = REPO / "test_data/hg_test.lengths"
bumbl = REPO / "test_data/hg_test.bumbl"
bi = Path("/tmp/hg_probe.bi")

if not bi.exists():
    subprocess.run(["shredtools", "index", str(bumbl), "--multi", "-o", str(bi)], check=True)

names = m.get_contig_names(str(lens))
lengths = m.get_sequence_lengths(str(lens), multilengths=True)
idx = sutils.parse_index(str(bi), seq_idx=0)

mums = m.MUMdata(str(bumbl), sort=False)
col = mums.starts[:, 0]
present = col >= 0
max_end = int((col[present] + mums.lengths[present]).max())

CANDIDATES = [
    ("chrM_empty", "CHM13#0#chrM:100-5000"),
    ("chr1_start", "CHM13#0#chr1:1-50000"),
    ("chr1_sparse_tail", "CHM13#0#chr1:248384423-248386328"),
    ("chrY_sparse_tail", "CHM13#0#chrY:62120000-62460000"),
    ("chr2_sparse_tail", "CHM13#0#chr2:242683100-242696700"),
    ("chrY_past_last_mum", "CHM13#0#chrY:62400000-62459000"),
]

print(f"max_mum_end={max_end:,} genome_total={sum(lengths[0]):,}\n")
for label, region in CANDIDATES:
    coords = sutils.convert_local_to_global_coords(region, names[0], lengths[0])
    gr = sutils.get_mum_ranges_region(idx, coords)
    gf = get_mum_ranges_flanks(idx, coords)
    print(f"{label}: {region}")
    print(f"  global={coords} region={gr is not None} flanks={gf is not None}")
    if gr is not None:
        print(f"  region_ranges={gr}")
    if gf is not None:
        print(f"  flank_ranges shape={gf.shape}")
