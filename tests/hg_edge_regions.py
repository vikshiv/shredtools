"""
Human (hg_test) edge regions derived from test_data/hg_test.lengths + bumbl.

Re-discover with: conda run -n shredtools python tests/_probe_hg_edges.py

Seq 0 = CHM13 (100000_CHM13.pri.fa), 25 contigs, total 3,117,292,070 bp.
Last MUM end on seq 0: 3,116,927,299 (364,771 bp before genome end; chrM has 0 MUMs).

| Region constant | Coordinates | Behavior |
|-----------------|-------------|----------|
| CHRM_EMPTY | chrM:100-5000 | No index bins (global past max_bin); subset empty exit 0; extract fail |
| CHRY_TAIL | chrY:62400000-62459000 | Past last MUM on chrY; no index bins |
| CHR1_START | chr1:1-50000 | Index bins exist; subset keeps 6 rows; extract flanks fail |
| CHR1_SPARSE_TAIL | chr1:248384423-248386328 | Index loads candidates; subset keeps 0; extract exit 0 but 0-line BED |
| CHR1_MID | chr1:1000000-2000000 | Normal query (sanity check) |
"""

# Empty contig — 0 MUMs anywhere on CHM13#0#chrM (len 16,569)
CHRM_EMPTY = "CHM13#0#chrM:100-5000"

# chrY tail: last MUM ends ~62,111,827; query entirely in 348 kb unmum gap
CHRY_TAIL = "CHM13#0#chrY:62400000-62459000"

# chr1 telomere: first MUM starts at global 1,704; flanks cannot bound query start
CHR1_START = "CHM13#0#chr1:1-50000"

# chr1 tail: last MUM local end 248,383,423; query in 3.9 kb gap before contig end
CHR1_SPARSE_TAIL = "CHM13#0#chr1:248384423-248386328"

# Well-covered interior region
CHR1_MID = "CHM13#0#chr1:1000000-2000000"

# chr2 tail with sparse gap but flanking bins still work for extract
CHR2_SPARSE_TAIL = "CHM13#0#chr2:242683100-242696700"
