#!/usr/bin/env python3
"""Synteny plot for shredtools-enhanced MUMs (globals + purple partial gap-fills).

Fork of mumemto untangle-branch viz_mums with:
  - global multi-MUMs as blue collinear blocks
  - partial MUMs (starts == -1) as purple collinear blocks on their present subset
  - optional genome y-order by partial co-occurrence (Fiedler seriation; --reorder)
  - retained --untangle / --mode spaced contig layout

Requires contig-aware mumemto.utils from the untangle tree for --untangle /
contig-coord inputs. Set MUMEMTO_UNTANGLE_ROOT to the package parent that
contains mumemto/utils.py (default: /homes/shivakumar/hps-vshiv/mumemto/mumemto).
"""

from __future__ import annotations

from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
import importlib.util
import os
import sys
import numpy as np
import argparse

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_UNTANGLE_ROOT = "/homes/shivakumar/hps-vshiv/mumemto/mumemto"


def _load_contig_utils():
    """Load contig-aware utils without importing mumemto package __init__ (C ext)."""
    roots = []
    env = os.environ.get("MUMEMTO_UNTANGLE_ROOT")
    if env:
        roots.append(env)
    roots.append(_DEFAULT_UNTANGLE_ROOT)

    for root in roots:
        utils_path = os.path.join(root, "mumemto", "utils.py")
        if not os.path.isfile(utils_path):
            continue
        spec = importlib.util.spec_from_file_location("_mumemto_utils_untangle", utils_path)
        mod = importlib.util.module_from_spec(spec)
        # archive/untangle.py does `from utils import ...`
        sys.modules["utils"] = mod
        sys.modules["_mumemto_utils_untangle"] = mod
        spec.loader.exec_module(mod)
        if not hasattr(mod.MUMdata, "has_contig_coords"):
            continue
        return mod

    # Fall back to installed mumemto (may lack contig helpers)
    try:
        from mumemto import utils as mod  # noqa: WPS433
        return mod
    except ImportError as e:
        raise ImportError(
            "Could not load mumemto.utils. Install mumemto or set "
            "MUMEMTO_UNTANGLE_ROOT to the untangle-branch package root "
            f"(tried: {roots})."
        ) from e


_utils = _load_contig_utils()
MUMdata = _utils.MUMdata
find_coll_blocks = _utils.find_coll_blocks
get_sequence_lengths = _utils.get_sequence_lengths
get_seq_paths = _utils.get_seq_paths
ensure_contig_coords = getattr(_utils, "ensure_contig_coords", None)

if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)
from untangle import (  # noqa: E402
    ContigLayout,
    untangle_layout,
    prepare_contig_plot,
    gapped_genome_length,
    apply_contig_layout,
    contig_local_to_concat,
)


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description=(
            "Synteny plot of enhanced MUMs: global multi-MUMs (blue) plus "
            "partial gap-fill MUMs as purple collinear blocks."
        )
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--input-prefix", "-i", dest="prefix",
        help="prefix for mums and lengths (uses *.mums over *.bumbl if both exist)",
    )
    group.add_argument(
        "--mums", "-m", dest="mumfile",
        help="path to *.mums (or *.bumbl) file",
    )

    parser.add_argument(
        "--lengths", "-l", dest="lens",
        help="lengths file (default: <prefix>.lengths)",
    )
    parser.add_argument(
        "--filelist", "-f", dest="filelist",
        help="filelist defines genome y-order (overrides --reorder). "
             "First col = full path from lengths file",
    )
    parser.add_argument(
        "--reorder", dest="reorder", action="store_true", default=False,
        help="auto-order genomes by partial-MUM co-occurrence (Fiedler seriation); "
             "off by default. Ignored if --filelist is given",
    )
    parser.add_argument(
        "--labels", "-n", dest="labels",
        help="labels for each sequence, one per line, matching filelist or lengths order",
    )
    parser.add_argument(
        "--len-filter", "-L", dest="lenfilter",
        help="only plot MUMs longer than threshold", default=0, type=int,
    )
    parser.add_argument(
        "--subsample", "-s", dest="subsample",
        help="subsample every Nth mum", default=1, type=int,
    )
    parser.add_argument(
        "--center", "-c", dest="center", action="store_true", default=False,
        help="center plot",
    )
    parser.add_argument(
        "--inversion-color", "-ic", dest="inv_color",
        help="color for inversions", default="green",
    )
    parser.add_argument(
        "--mum-color", "-mc", dest="mum_color",
        help="color for global multi-MUMs (default: #00A2FF)",
        default="#00A2FF", type=str,
    )
    parser.add_argument(
        "--partial-color", "-pc", dest="partial_color",
        help="color for partial MUM blocks (default: #9B30FF)",
        default="#9B30FF", type=str,
    )
    parser.add_argument(
        "--alpha", "-a", dest="alpha",
        help="opacity of mums [0-1] (default: 0.8 for blocks, 0.05 for MUMs)",
        type=float,
    )
    parser.add_argument(
        "--linewidth", "-lw", dest="linewidth",
        help="linewidth of mums (default: 0 for blocks, 0.05 for MUMs)",
        type=float,
    )
    parser.add_argument(
        "--fout", "-o", dest="filename",
        help="plot fname (default: input_prefix)",
    )
    parser.add_argument(
        "--dims", dest="size",
        help="fig dimensions (inches) (default: 6.4, 4.8)",
        default=(6.4, 4.8), type=float, nargs=2,
    )
    parser.add_argument(
        "--region", dest="region",
        help="genomic region to zoom in (global/concat coords)",
        default=None, type=float, nargs=2,
    )
    parser.add_argument("--dpi", "-d", dest="dpi", help="dpi", default=500, type=int)
    parser.add_argument(
        "--verbose", "-v", dest="verbose", action="store_true", default=False,
    )
    parser.add_argument(
        "--no-coll-block", "-b", dest="no_coll_block",
        help="plot individual MUMs instead of collinear blocks",
        action="store_true", default=False,
    )
    parser.add_argument(
        "--max-gap-len", "-g", dest="max_break",
        help="max break between collinear mums within a block (default: <1px)",
        default=None, type=int,
    )

    parser.add_argument(
        "--mode", dest="mode",
        choices=["normal", "delineated", "gapped", "spaced"],
        help=(
            "Plotting mode for multi-fasta inputs: normal=concat; "
            "delineated=concat with red contig ticks; "
            "gapped=column-aligned gaps (equal contig counts); "
            "spaced=per-genome contig gaps + scale bar (default with --untangle)"
        ),
        default="normal",
    )
    parser.add_argument(
        "--spacer", dest="spacer", type=float, default=0.1,
        help=(
            "Inter-contig gap for --mode spaced/gapped, as a fraction of the largest "
            "contig length (default: 0.1)."
        ),
    )
    parser.add_argument(
        "--untangle", dest="untangle", action="store_true", default=False,
        help=(
            "Reorder/orient contigs vs genome 0 (keeps ref fixed; others by median "
            "mapped position / majority strand). Implies --mode spaced. "
            "Requires contig coords or multilengths."
        ),
    )

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    if args.mumfile:
        base, ext = os.path.splitext(args.mumfile)
        if ext in (".mums", ".bumbl"):
            args.prefix = base
        else:
            args.prefix = args.mumfile
    elif args.prefix:
        if args.prefix.endswith(".mums") or args.prefix.endswith(".bumbl"):
            args.prefix = os.path.splitext(args.prefix)[0]
        if os.path.exists(args.prefix + ".mums"):
            args.mumfile = args.prefix + ".mums"
        elif os.path.exists(args.prefix + ".bumbl"):
            args.mumfile = args.prefix + ".bumbl"
        else:
            parser.error("No *.mums or *.bumbl file found for prefix %s" % args.prefix)
    else:
        parser.error("Either --mums or --prefix must be provided")

    if args.lens is None:
        args.lens = args.prefix + ".lengths"
        if not os.path.exists(args.lens):
            raise FileNotFoundError(
                f"Lengths file {args.lens} not found, and no lengths file provided"
            )

    if not args.filename:
        args.filename = args.prefix

    if args.alpha is None:
        args.alpha = 0.05 if args.no_coll_block else 0.8
    if args.linewidth is None:
        args.linewidth = 0.05 if args.no_coll_block else 0

    if args.untangle and args.mode == "normal":
        args.mode = "spaced"

    return args


def _nice_scale_bp(visible_span_bp):
    span = float(visible_span_bp)
    if span <= 0:
        return 1, "1 bp"
    candidates = [
        (1, "1 bp"), (2, "2 bp"), (5, "5 bp"),
        (10, "10 bp"), (20, "20 bp"), (50, "50 bp"),
        (100, "100 bp"), (200, "200 bp"), (500, "500 bp"),
        (1e3, "1 kb"), (2e3, "2 kb"), (5e3, "5 kb"),
        (1e4, "10 kb"), (2e4, "20 kb"), (5e4, "50 kb"),
        (1e5, "100 kb"), (2e5, "200 kb"), (5e5, "500 kb"),
        (1e6, "1 Mb"), (2e6, "2 Mb"), (5e6, "5 Mb"),
        (1e7, "10 Mb"), (2e7, "20 Mb"), (5e7, "50 Mb"),
        (1e8, "100 Mb"), (2e8, "200 Mb"), (5e8, "500 Mb"),
        (1e9, "1 Gb"), (2e9, "2 Gb"), (5e9, "5 Gb"),
    ]
    target = span * 0.12
    fits = [c for c in candidates if 0.05 * span <= c[0] <= 0.22 * span]
    pool = fits if fits else candidates
    best = min(pool, key=lambda c: abs(c[0] - target))
    return int(best[0]), best[1]


def add_scale_bar(ax, color="black"):
    x0, x1 = ax.get_xlim()
    span = abs(x1 - x0)
    if span <= 0:
        return
    length_bp, label = _nice_scale_bp(span)
    bar_frac = length_bp / span
    right = 0.98
    left = max(0.02, right - bar_frac)
    y = -0.12
    tick = 0.015
    ax.plot(
        [left, right], [y, y], transform=ax.transAxes,
        color=color, linewidth=1.5, solid_capstyle="butt", clip_on=False,
    )
    ax.plot(
        [left, left], [y - tick, y + tick], transform=ax.transAxes,
        color=color, linewidth=1.2, clip_on=False,
    )
    ax.plot(
        [right, right], [y - tick, y + tick], transform=ax.transAxes,
        color=color, linewidth=1.2, clip_on=False,
    )
    ax.text(
        (left + right) / 2, y - 0.02, label, transform=ax.transAxes,
        ha="center", va="top", fontsize=8, color=color, clip_on=False,
    )


def points_to_poly(points):
    starts, ends = tuple(zip(*points))
    points = starts + ends[::-1]
    return points


def get_mum_polygons(mums, centering, color="#00A2FF", inv_color="red"):
    polygons = []
    colors = []
    for (l, starts, strands) in mums:
        inverted = not strands[0]
        points = []
        for idx, (x, strand) in enumerate(zip(starts, strands)):
            if x == -1:
                if len(points) > 2:
                    polygons.append(points_to_poly(points[:-1]))
                    colors.append(color)
                continue
            points.append(((centering[idx] + x, idx), (centering[idx] + x + l, idx)))
            if not inverted and not strand:
                inverted = True
                if len(points) > 2:
                    polygons.append(points_to_poly(points[:-1]))
                    colors.append(color)
                polygons.append(points_to_poly(points[-2:]))
                colors.append(inv_color)
                points = [points[-1]]
            elif inverted and strand:
                inverted = False
                if len(points) > 2:
                    polygons.append(points_to_poly(points[:-1]))
                    colors.append(color)
                polygons.append(points_to_poly(points[-2:]))
                colors.append(inv_color)
                points = [points[-1]]
        if len(points) >= 2:
            polygons.append(points_to_poly(points))
            colors.append(color)
    return polygons, colors


def get_block_polygons(collinear_blocks, mums, centering, color="#00A2FF", inv_color="red"):
    polygons = []
    colors = []
    for _, (l, r) in enumerate(collinear_blocks):
        strands = mums[l].strands
        inverted = not strands[0]
        points = []
        left, right = mums[l].starts.copy(), (mums[r].starts + mums[r].length).copy()
        for idx in range(len(left)):
            if left[idx] != -1 and right[idx] != -1 and left[idx] > right[idx]:
                left[idx], right[idx] = right[idx], left[idx]
        for idx, strand in enumerate(strands):
            if left[idx] == -1:
                if len(points) > 2:
                    polygons.append(points_to_poly(points))
                    colors.append(color)
                points = []
                continue
            points.append(
                ((centering[idx] + left[idx], idx), (centering[idx] + right[idx], idx))
            )
            if not inverted and not strand:
                inverted = True
                if len(points) > 2:
                    polygons.append(points_to_poly(points[:-1]))
                    colors.append(color)
                polygons.append(points_to_poly(points[-2:]))
                colors.append(inv_color)
                points = [points[-1]]
            elif inverted and strand:
                inverted = False
                if len(points) > 2:
                    polygons.append(points_to_poly(points[:-1]))
                    colors.append(color)
                polygons.append(points_to_poly(points[-2:]))
                colors.append(inv_color)
                points = [points[-1]]
        if len(points) >= 2:
            polygons.append(points_to_poly(points))
            colors.append(color)
    return polygons, colors


def plot(args, genome_lengths, polygons, colors, centering, dpi=500, size=None, genomes=None, filename=None):
    fig, ax = plt.subplots()
    max_length = max(genome_lengths) if genome_lengths else 0
    spacer_bp = int(getattr(args, "spacer_bp", 0) or 0)

    if args.mode == "normal" or args.mode == "delineated":
        for idx, g in enumerate(genome_lengths):
            ax.plot(
                [centering[idx] + 0, centering[idx] + g], [idx, idx],
                alpha=0.2, linewidth=0.75, c="black",
            )
    elif args.mode == "spaced":
        for idx, ml in enumerate(args.multilengths):
            pos = centering[idx]
            for i, L in enumerate(ml):
                L = int(L)
                if L > 0:
                    ax.plot(
                        [pos, pos + L], [idx, idx],
                        alpha=0.35, linewidth=1.0, c="black", solid_capstyle="butt",
                    )
                pos += L
                if i < len(ml) - 1:
                    pos += spacer_bp
    elif args.mode == "gapped":
        offsets = np.array(
            [0] + (args.multilengths.max(axis=0) + args.spacer).cumsum().tolist()[:-1]
        )
        vert_seps = [p - (args.spacer) / 2 for p in offsets] + [
            args.multilengths.max(axis=0).sum()
            + args.spacer * (args.multilengths.shape[1] - 1)
        ]
        for p in vert_seps[1:-1]:
            ax.plot([p, p], [0, len(genome_lengths) - 1], alpha=0.5, linewidth=1, color="black")
        for idx in range(args.multilengths.shape[0]):
            for i, offset in enumerate(args.multilengths[idx]):
                ax.plot(
                    [centering[idx] + offsets[i], centering[idx] + offsets[i] + offset],
                    [idx, idx], alpha=0.2, linewidth=0.25,
                )
        chr_markers = [
            vert_seps[idx - 1] + ((vert_seps[idx] - vert_seps[idx - 1]) / 2)
            for idx in range(1, len(vert_seps))
        ]
        ax.set_xticks(chr_markers)
        ax.set_xticklabels(range(1, len(chr_markers) + 1))
    else:
        for idx, g in enumerate(genome_lengths):
            ax.plot(
                [centering[idx] + 0, centering[idx] + g], [idx, idx],
                alpha=0.2, linewidth=0.75,
            )

    if polygons:
        ax.add_collection(
            PolyCollection(
                polygons, linewidths=args.linewidth, alpha=args.alpha,
                edgecolors=colors, facecolors=colors,
            )
        )

    if args.mode == "delineated":
        for idx in range(len(args.multilengths)):
            cur_offsets = np.cumsum(args.multilengths[idx])
            for i, offset in enumerate([0] + cur_offsets[:-1].tolist()):
                ax.plot(
                    [centering[idx] + offset, centering[idx] + offset],
                    [idx - 0.25, idx + 0.25], alpha=1, linewidth=0.25, color="red",
                )

    ax.yaxis.set_ticks(list(range(len(genome_lengths))))
    ax.tick_params(axis="y", which="both", length=0)
    if genomes:
        ax.set_yticklabels(genomes, fontsize=8)
    else:
        ax.yaxis.set_ticklabels([])
    if args.mode == "gapped":
        ax.set_xlabel("chromosome")
    elif args.mode == "spaced":
        ax.set_xlabel("")
        ax.set_xticks([])
        ax.tick_params(axis="x", which="both", length=0, labelbottom=False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    else:
        ax.set_xlabel("genomic position")
    ax.set_ylabel("sequences")
    ax.set_ylim(-0.25, len(genome_lengths) - 1 + 0.25)
    if args.mode == "gapped":
        ax.set_xlim(
            0,
            args.multilengths.max(axis=0).sum()
            + args.spacer * (args.multilengths.shape[1] - 1),
        )
    elif args.region is not None:
        ax.set_xlim(*args.region)
    else:
        ax.set_xlim(0, max_length)
    ax.invert_yaxis()
    if args.mode == "spaced":
        add_scale_bar(ax)
        fig.subplots_adjust(bottom=0.18)
    else:
        fig.set_tight_layout(True)

    if size:
        fig.set_size_inches(*size)
    if filename:
        if os.path.splitext(filename)[1] not in [".png", ".pdf", ".svg"]:
            filename = filename + ".png"
        if os.path.dirname(filename):
            path = filename
        else:
            path = os.path.join(os.path.dirname(args.mumfile), filename)
        fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.25)
    return ax


def offset_mums(args, mums, spacer=None, blocks=None):
    offset = args.multilengths
    if spacer is None:
        spacer = args.spacer
    NUM_SEQS = len(offset)
    offsets = np.cumsum(offset, axis=1)
    contig_idx = np.array(
        [
            np.searchsorted(offsets[idx], mums.starts[:, idx], side="right")
            for idx in range(NUM_SEQS)
        ]
    ).transpose()
    if blocks is not None:
        new_blocks = []
        for l, r in blocks:
            if not np.all(contig_idx[l, :] == contig_idx[r, :]):
                if l == r - 1:
                    continue
                diffs = ~(np.diff(contig_idx[l : r + 1, :], axis=0) == 0).all(axis=1)
                old_l = l
                for s in np.where(diffs)[0] + 1:
                    if s == 1:
                        l = old_l + 2
                        continue
                    cur = old_l + s
                    if l < cur - 1 and cur - 1 <= r:
                        new_blocks.append((l, cur - 1))
                        l = cur + 1
                if l < r:
                    new_blocks.append((l, r))
            else:
                new_blocks.append((l, r))
    else:
        new_blocks = None
    left_start = np.hstack((np.zeros((offsets.shape[0], 1), dtype=int), offsets[:, :-1]))
    rel_offsets = mums.starts - left_start[np.arange(NUM_SEQS), contig_idx]
    partial_mask = mums.starts != -1
    new_starts = (
        np.array([0] + (offset.max(axis=0) + spacer).cumsum().tolist()[:-1])[contig_idx]
        + rel_offsets
    )
    mums.starts[partial_mask] = new_starts[partial_mask]
    return new_blocks


def insert_spaced_spacers(mums, multilengths, spacer, blocks=None):
    spacer = int(spacer)
    n_seqs = len(multilengths)
    contig_idx = np.zeros((mums.num_mums, n_seqs), dtype=np.int32)
    for g, ml in enumerate(multilengths):
        offsets = np.cumsum(ml)
        starts = mums.starts[:, g]
        present = starts >= 0
        cid = np.zeros(mums.num_mums, dtype=np.int32)
        cid[present] = np.searchsorted(offsets, starts[present], side="right")
        cid[present] = np.minimum(cid[present], len(ml) - 1)
        contig_idx[:, g] = cid
        if spacer:
            mums.starts[present, g] = starts[present] + cid[present] * spacer

    if blocks is None:
        return None
    new_blocks = []
    for l, r in blocks:
        if not np.all(contig_idx[l, :] == contig_idx[r, :]):
            if l == r - 1:
                continue
            diffs = ~(np.diff(contig_idx[l : r + 1, :], axis=0) == 0).all(axis=1)
            old_l = l
            cur_l = l
            for s in np.where(diffs)[0] + 1:
                if s == 1:
                    cur_l = old_l + 2
                    continue
                cur = old_l + s
                if cur_l < cur - 1 and cur - 1 <= r:
                    new_blocks.append((cur_l, cur - 1))
                    cur_l = cur + 1
            if cur_l < r:
                new_blocks.append((cur_l, r))
        else:
            new_blocks.append((l, r))
    return new_blocks


# ---- Enhanced partial-MUM helpers ------------------------------------------

def presence_mask_key(starts_row):
    return tuple(bool(x != -1) for x in starts_row)


def iter_partial_stretches(starts):
    """Yield (start_idx, end_idx_exclusive, presence_tuple) for contiguous partial runs."""
    n = starts.shape[0]
    i = 0
    while i < n:
        if not np.any(starts[i] == -1):
            i += 1
            continue
        mask = presence_mask_key(starts[i])
        j = i + 1
        while j < n and np.any(starts[j] == -1) and presence_mask_key(starts[j]) == mask:
            j += 1
        yield i, j, mask
        i = j


def assign_stretch_ids(starts):
    """Integer stretch id per MUM row; -1 for strict MUMs."""
    ids = np.full(starts.shape[0], -1, dtype=np.int32)
    sid = 0
    for i, j, _mask in iter_partial_stretches(starts):
        ids[i:j] = sid
        sid += 1
    return ids


def fiedler_genome_order(starts, lengths, n_seqs):
    """Order genomes by Fiedler vector of length-weighted partial co-occurrence."""
    affinity = np.zeros((n_seqs, n_seqs), dtype=np.float64)
    for i, j, mask in iter_partial_stretches(starts):
        present = [g for g, p in enumerate(mask) if p]
        if len(present) < 2:
            continue
        w = float(np.sum(lengths[i:j]))
        for a_i, a in enumerate(present):
            for b in present[a_i + 1 :]:
                affinity[a, b] += w
                affinity[b, a] += w

    if not np.any(affinity):
        return np.arange(n_seqs, dtype=int)

    degree = affinity.sum(axis=1)
    # Isolated genomes keep identity contribution so Laplacian is well-defined
    deg_safe = np.where(degree > 0, degree, 1.0)
    d_inv_sqrt = 1.0 / np.sqrt(deg_safe)
    # Normalized Laplacian L = I - D^{-1/2} A D^{-1/2}
    scaled = affinity * d_inv_sqrt[:, None] * d_inv_sqrt[None, :]
    lap = np.eye(n_seqs) - scaled
    # Zero out rows/cols for genomes never in a partial stretch (keep near end/start stably)
    unused = degree == 0
    if np.all(unused):
        return np.arange(n_seqs, dtype=int)
    used = np.where(~unused)[0]
    if len(used) <= 2:
        return np.arange(n_seqs, dtype=int)

    sub = lap[np.ix_(used, used)]
    evals, evecs = np.linalg.eigh(sub)
    # Second-smallest eigenvector
    fiedler = evecs[:, 1] if sub.shape[0] > 1 else evecs[:, 0]
    order_used = used[np.argsort(fiedler)]
    order = list(order_used)
    # Append unused genomes in original relative order
    for g in range(n_seqs):
        if unused[g]:
            order.append(g)
    return np.asarray(order, dtype=int)


def apply_genome_order(mums, order, seq_lengths, multilengths, mode):
    order = np.asarray(order, dtype=int)
    mums.starts = mums.starts[:, order]
    mums.strands = mums.strands[:, order]
    if getattr(mums, "contig_ids", None) is not None:
        mums.contig_ids = mums.contig_ids[:, order]
    seq_lengths = [seq_lengths[i] for i in order]
    if multilengths is not None:
        if mode == "gapped" and isinstance(multilengths, np.ndarray):
            multilengths = multilengths[order]
        else:
            multilengths = [multilengths[i] for i in order]
    return mums, seq_lengths, multilengths


def partial_block_polygons(partials, stretch_ids, centering, max_break, color, inv_color, verbose=False):
    """Find collinear blocks within each stretch (projected to present genomes) and polygonize."""
    polygons, colors = [], []
    if partials is None or partials.num_mums == 0:
        return polygons, colors

    unique_ids = [sid for sid in np.unique(stretch_ids) if sid >= 0]
    n_blocks = 0
    for sid in unique_ids:
        rows = np.where(stretch_ids == sid)[0]
        stretch = partials[rows]
        present = np.where(np.any(stretch.starts != -1, axis=0))[0]
        if len(present) < 2 or stretch.num_mums == 0:
            continue
        projected = stretch[:, present]
        # Collinear on the present subset only
        blocks = find_coll_blocks(projected, max_break=max_break, verbose=False)
        if not blocks:
            continue
        n_blocks += len(blocks)
        # Map block indices (into stretch) — projected shares row order with stretch
        polys, cols = get_block_polygons(
            blocks, stretch, centering, color=color, inv_color=inv_color,
        )
        polygons.extend(polys)
        colors.extend(cols)

    if verbose:
        print(
            f"partial stretches: {len(unique_ids)}; purple collinear blocks: {n_blocks}",
            file=sys.stderr,
        )
    return polygons, colors


def apply_contig_layout_to_partials(partials, multilengths, layout, spacer=0):
    """Apply the same ContigLayout used for strict MUMs to partials; return concat coords."""
    if partials is None or partials.num_mums == 0:
        return partials, multilengths
    try:
        laid, ml_new = apply_contig_layout(partials, multilengths, layout)
        plot_mums = contig_local_to_concat(
            laid, ml_new, keep_contig_ids=False, spacer=int(spacer),
        )
        return plot_mums, ml_new
    except Exception as exc:  # noqa: BLE001 — degrade gracefully per plan
        print(
            f"Warning: could not apply contig layout to partials ({exc}); "
            "plotting partials without flips.",
            file=sys.stderr,
        )
        try:
            plot_mums = contig_local_to_concat(
                partials, multilengths, keep_contig_ids=False, spacer=int(spacer),
            )
            return plot_mums, multilengths
        except Exception as exc2:  # noqa: BLE001
            print(
                f"Warning: partial contig concat failed ({exc2}); skipping purple layer.",
                file=sys.stderr,
            )
            return None, multilengths


def default_max_break(seq_lengths, dpi, width_in):
    max_length = max(seq_lengths) if seq_lengths else 1
    bp_per_inch = max_length / (dpi * width_in)
    return int(min(bp_per_inch, 100000))


def main(args):
    if args.mode != "normal" or args.untangle:
        try:
            offset = get_sequence_lengths(args.lens, multilengths=True)
            seq_lengths = [sum(o) for o in offset]
            if args.mode == "gapped" and len(set([len(o) for o in offset])) > 1:
                print(
                    "Warning: gapped mode requires the same number of sequences per "
                    "input FASTA file. Using spaced mode instead.",
                    file=sys.stderr,
                )
                args.mode = "spaced"
                args.multilengths = offset
            elif args.mode == "gapped":
                args.multilengths = np.array(offset)
            else:
                args.multilengths = offset
        except ValueError:
            if args.untangle:
                print(
                    "Error: --untangle requires multilengths in %s" % args.lens,
                    file=sys.stderr,
                )
                sys.exit(1)
            print(
                "Warning: Multi-FASTA lengths not available in %s. Treating input "
                "FASTAs as a single sequence instead." % args.lens,
                file=sys.stderr,
            )
            args.mode = "normal"
            seq_lengths = get_sequence_lengths(args.lens)
            args.multilengths = None
    else:
        seq_lengths = get_sequence_lengths(args.lens)
        args.multilengths = None

    args.spacer_bp = 0
    if args.mode == "gapped":
        args.spacer = args.spacer * args.multilengths.max(axis=0).max()
        args.spacer_bp = int(args.spacer)
    elif args.mode == "spaced" and args.multilengths is not None:
        max_contig = max(max(ml) for ml in args.multilengths if len(ml))
        args.spacer_bp = int(args.spacer * max_contig)

    mums = MUMdata(
        args.mumfile,
        lenfilter=args.lenfilter,
        subsample=args.subsample,
        verbose=args.verbose,
    )
    if args.verbose:
        has_cc = getattr(mums, "has_contig_coords", False)
        if callable(has_cc):
            has_cc = mums.has_contig_coords
        print(
            f"Found {mums.num_mums} MUMs" + (" (contig coords)" if has_cc else ""),
            file=sys.stderr,
        )

    n_seqs = mums.num_seqs
    genome_order = None

    if args.filelist:
        cur_order = get_seq_paths(args.lens)
        file_order = [l.split()[0] for l in open(args.filelist, "r").read().splitlines()]
        try:
            genome_order = np.array([cur_order.index(l) for l in file_order])
        except ValueError:
            print(
                "Error: sequence in filelist not found in lengths file. Ensure the "
                "filelist paths match the lengths file paths exactly.",
                file=sys.stderr,
            )
            sys.exit(1)
    elif args.reorder:
        genome_order = fiedler_genome_order(mums.starts, mums.lengths, n_seqs)
        if args.verbose and not np.array_equal(genome_order, np.arange(n_seqs)):
            print(f"auto genome order: {genome_order.tolist()}", file=sys.stderr)

    if genome_order is not None:
        mums, seq_lengths, args.multilengths = apply_genome_order(
            mums, genome_order, seq_lengths, args.multilengths, args.mode,
        )

    if args.labels:
        if args.labels.endswith(".lengths"):
            genome_names = [
                os.path.splitext(os.path.basename(l))[0]
                for l in get_seq_paths(args.labels)
            ]
            if genome_order is not None:
                genome_names = [genome_names[i] for i in genome_order]
        else:
            genome_names = open(args.labels, "r").read().splitlines()
    else:
        genome_names = None

    # Split strict vs partial before any filter_pmums
    is_partial = np.any(mums.starts == -1, axis=1)
    stretch_ids_full = assign_stretch_ids(mums.starts)
    if np.any(~is_partial):
        strict = mums[~is_partial]
    else:
        cids = None
        if getattr(mums, "contig_ids", None) is not None:
            cids = np.zeros((0, n_seqs), dtype=np.int32)
        strict = MUMdata.from_arrays(
            np.array([], dtype=np.uint32),
            np.zeros((0, n_seqs), dtype=np.int64),
            np.zeros((0, n_seqs), dtype=bool),
            contig_ids=cids,
        )
    partials = mums[is_partial] if np.any(is_partial) else None
    partial_stretch_ids = stretch_ids_full[is_partial] if np.any(is_partial) else None

    if args.verbose:
        n_partial = int(np.sum(is_partial))
        print(
            f"strict MUMs: {mums.num_mums - n_partial}; partial MUMs: {n_partial}",
            file=sys.stderr,
        )

    has_contig = getattr(mums, "has_contig_coords", False)
    if callable(has_contig):
        has_contig = mums.has_contig_coords
    contig_mode = bool(args.untangle or has_contig)

    if contig_mode:
        if ensure_contig_coords is None:
            print(
                "Error: contig coordinates / --untangle require contig-aware "
                "mumemto.utils (set MUMEMTO_UNTANGLE_ROOT).",
                file=sys.stderr,
            )
            sys.exit(1)
        if args.multilengths is None:
            print(
                "Error: contig coordinates / --untangle require multilengths",
                file=sys.stderr,
            )
            sys.exit(1)

        ml_orig = [list(ml) for ml in args.multilengths]
        ensure_contig_coords(strict, ml_orig)
        if partials is not None and partials.num_mums > 0:
            ensure_contig_coords(partials, ml_orig)

        if strict.num_mums == 0 and (partials is None or partials.num_mums == 0):
            print("No MUMs found after split.", file=sys.stderr)
            return

        layout = ContigLayout.identity(ml_orig)
        if args.untangle and strict.num_mums > 0:
            layout = untangle_layout(strict, ml_orig, ref=0, verbose=args.verbose)

        if args.max_break is None:
            args.max_break = default_max_break(seq_lengths, args.dpi, args.size[0])
        if args.mode == "spaced":
            max_contig = max(max(ml) for ml in ml_orig if len(ml))
            args.spacer_bp = int(args.spacer * max_contig)

        spacer = args.spacer_bp if args.mode == "spaced" else 0
        if strict.num_mums > 0:
            strict, blocks, args.multilengths = prepare_contig_plot(
                strict, ml_orig, layout,
                max_break=args.max_break, verbose=args.verbose, spacer=spacer,
            )
            strict.blocks = blocks
        else:
            args.multilengths = layout.multilengths_reordered(ml_orig)

        if partials is not None and partials.num_mums > 0:
            partials, _ml = apply_contig_layout_to_partials(
                partials, ml_orig, layout, spacer=spacer,
            )
            if partials is None:
                partial_stretch_ids = None

        if args.mode == "spaced":
            seq_lengths = [
                gapped_genome_length(ml, args.spacer_bp) for ml in args.multilengths
            ]
        else:
            seq_lengths = [sum(ml) for ml in args.multilengths]
        if args.mode == "normal":
            args.mode = "spaced"
        if args.verbose and strict.num_mums > 0:
            nblk = len(getattr(strict, "blocks", None) or [])
            print(
                f"contig layout applied; {strict.num_mums} strict MUMs, {nblk} blocks"
                + (f"; spacer={args.spacer_bp}bp" if args.mode == "spaced" else ""),
                file=sys.stderr,
            )

    if args.mode == "spaced" and not contig_mode and args.multilengths is not None:
        seq_lengths = [
            gapped_genome_length(ml, args.spacer_bp) for ml in args.multilengths
        ]

    max_length = max(seq_lengths) if seq_lengths else 0
    centering = [0] * len(seq_lengths)
    if args.center:
        centering = [(max_length - g) / 2 for g in seq_lengths]

    if args.max_break is None:
        args.max_break = default_max_break(seq_lengths, args.dpi, args.size[0])

    polygons, colors = [], []

    if args.no_coll_block:
        if not contig_mode:
            if args.mode == "gapped":
                if strict.num_mums > 0:
                    offset_mums(args, strict)
                if partials is not None and partials.num_mums > 0:
                    offset_mums(args, partials)
            elif args.mode == "spaced" and args.multilengths is not None:
                if strict.num_mums > 0:
                    insert_spaced_spacers(strict, args.multilengths, args.spacer_bp)
                if partials is not None and partials.num_mums > 0:
                    insert_spaced_spacers(partials, args.multilengths, args.spacer_bp)

        if strict.num_mums > 0:
            p, c = get_mum_polygons(
                strict, centering, color=args.mum_color, inv_color=args.inv_color,
            )
            polygons.extend(p)
            colors.extend(c)
        if partials is not None and partials.num_mums > 0:
            p, c = get_mum_polygons(
                partials, centering, color=args.partial_color, inv_color=args.inv_color,
            )
            polygons.extend(p)
            colors.extend(c)
    else:
        if strict.num_mums > 0:
            if strict.blocks is None:
                collinear_blocks = find_coll_blocks(
                    strict, max_break=args.max_break, verbose=args.verbose,
                )
                if args.verbose:
                    print(
                        f"found {len(collinear_blocks)} global collinear blocks",
                        file=sys.stderr,
                    )
            else:
                collinear_blocks = strict.blocks
                if args.verbose:
                    print(
                        f"Using pre-computed collinear blocks: {len(collinear_blocks)} blocks",
                        file=sys.stderr,
                    )
            if not contig_mode:
                if args.mode == "gapped":
                    collinear_blocks = offset_mums(
                        args, strict, spacer=args.spacer, blocks=collinear_blocks,
                    )
                elif args.mode == "spaced" and args.multilengths is not None:
                    collinear_blocks = insert_spaced_spacers(
                        strict, args.multilengths, args.spacer_bp, blocks=collinear_blocks,
                    )
            p, c = get_block_polygons(
                collinear_blocks, strict, centering,
                color=args.mum_color, inv_color=args.inv_color,
            )
            polygons.extend(p)
            colors.extend(c)
        elif args.verbose:
            print("No strict MUMs to plot (blue layer empty).", file=sys.stderr)

        if partials is not None and partials.num_mums > 0 and partial_stretch_ids is not None:
            if not contig_mode:
                if args.mode == "gapped":
                    offset_mums(args, partials)
                elif args.mode == "spaced" and args.multilengths is not None:
                    insert_spaced_spacers(partials, args.multilengths, args.spacer_bp)
            p, c = partial_block_polygons(
                partials, partial_stretch_ids, centering, args.max_break,
                color=args.partial_color, inv_color=args.inv_color, verbose=args.verbose,
            )
            polygons.extend(p)
            colors.extend(c)

    if args.verbose:
        print("Rendering plot...", file=sys.stderr)
    plot(
        args, seq_lengths, polygons, colors, centering,
        genomes=genome_names, filename=args.filename, dpi=args.dpi, size=args.size,
    )
    if args.verbose:
        print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main(parse_arguments())
