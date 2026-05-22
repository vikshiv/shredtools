"""Extract-window MUM synteny plot (aligned with shredtools/extract_from_mums --plot)."""

from __future__ import annotations

import io
from typing import Sequence

import viz_mums


def plot(genome_lengths, polygons, colors, centering, xlims=None, size=None, genomes=None):
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.collections import PolyCollection
    from matplotlib import pyplot as plt

    fig, ax = plt.subplots()
    max_length = max(genome_lengths)
    for idx, g in enumerate(genome_lengths):
        ax.plot(
            [centering[idx] + 0, centering[idx] + g],
            [idx, idx],
            alpha=0.2,
            linewidth=0.75,
            c="black",
        )

    if xlims is not None:
        ax.set_xlim(*xlims)
    else:
        ax.set_xlim(0, max_length)
    ax.add_collection(
        PolyCollection(
            polygons,
            linewidths=0,
            alpha=0.8,
            edgecolors=colors,
            facecolors=colors,
        )
    )

    ax.yaxis.set_ticks(list(range(len(genome_lengths))))
    ax.tick_params(axis="y", which="both", length=0)
    if genomes:
        ax.set_yticklabels(genomes)
    else:
        ax.yaxis.set_ticklabels([])

    ax.set_xlabel("genomic position")
    ax.set_ylabel("sequences")
    ax.set_ylim(-0.25, len(genome_lengths) - 1 + 0.25)
    ax.invert_yaxis()
    fig.set_tight_layout(True)
    if size:
        fig.set_size_inches(*size)
    return fig, ax


def plot_extract(
    coords,
    mums,
    mum_bounds,
    other_coords,
    seq_idx,
    sequences: Sequence[int],
    seq_lengths: Sequence[int],
    genome_labels: Sequence[str] | None = None,
) -> bytes:
    """
    Render extract synteny to PNG bytes (window-relative coordinates).
    ``sequences`` lists source sequence indices; plot rows are 0..len(sequences)-1.
    """
    i0, i1 = int(mum_bounds[0]), int(mum_bounds[1])
    seq_list = [int(s) for s in sequences]
    plot_mums = mums.slice_rows(i0, i1, seq_list)
    for row_i, src_seq in enumerate(seq_list):
        plot_mums.offset_starts_col(row_i, int(other_coords[src_seq][0]))

    ref_row = seq_list.index(int(seq_idx))
    ref_offset = int(other_coords[seq_idx][0])
    start = int(coords[0]) - ref_offset
    end = int(coords[1]) - ref_offset

    centering = [0] * len(seq_list)
    poly, colors = viz_mums.get_mum_polygons(plot_mums, centering, inv_color="green")

    oc_sub = [other_coords[s] for s in seq_list]
    x_max = max(int(b) - int(a) for a, b in oc_sub)
    n_rows = len(seq_list)
    height = max(5.0, 0.12 * n_rows)
    genome_lengths = [int(seq_lengths[s]) for s in seq_list]
    labels = list(genome_labels) if genome_labels else None

    fig, ax = plot(
        genome_lengths,
        poly,
        colors,
        centering,
        xlims=(0, x_max),
        size=(10, height),
        genomes=labels,
    )
    ax.plot(
        [start, start],
        [ref_row - 0.5, ref_row + 0.5],
        color="red",
        linestyle="--",
        linewidth=1,
    )
    ax.plot(
        [end, end],
        [ref_row - 0.5, ref_row + 0.5],
        color="red",
        linestyle="--",
        linewidth=1,
    )

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    from matplotlib import pyplot as plt

    plt.close(fig)
    return buf.getvalue()
