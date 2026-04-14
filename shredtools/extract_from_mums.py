#!/usr/bin/env python3

import argparse
import sys, os
import numpy as np
from mumemto.utils import MUMdata, get_sequence_lengths, get_seq_paths, find_coll_blocks
import mumemto.viz_mums
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
from tqdm.auto import tqdm
from bisect import bisect_left

def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description="Shred a MUM file into smaller fragments and optionally visualize.")
    parser.add_argument("mum_file", type=str, help="Path to the input MUM file.")
    parser.add_argument("--seq-idx", '-s', type=int, help="Index of sequence with target region", required=True)
    parser.add_argument("--range", '-r', type=str, help="region coordinates. Format -> chr:start-end", required=True)
    parser.add_argument("--plot", action="store_true", help="If set, generate a visualization of the extracted region using MUMs from MUM file.")
    parser.add_argument("--plot-full", action="store_true", help="If set, generate a visualization of the extracted region using all MUMs from MUM file.")
    parser.add_argument("--fasta", action="store_true", help="If set, generate a FASTA file for each sequence that contains the target sequence. Otherwise, only write a BED file of coordinates.")
    parser.add_argument("--output", '-o', type=str, default="output", help="Output prefix for shreds. With --fasta, specifies directory to store shred sequences")
    parser.add_argument("--sequences", '-x', type=int, nargs='*', default=None, help="One or more sequence indices to output BED or FASTA for. By default, all sequences are included.")
    parser.add_argument('--lengths','-l', dest='lens', help='lengths file, first column is seq length in order of filelist')
    
    args = parser.parse_args(args)
    
    if not os.path.exists(args.mum_file):
        print(f"MUM file {args.mum_file} not found", file=sys.stderr)
        sys.exit(1)
    
    if args.lens is None:
        args.lens = os.path.splitext(args.mum_file)[0] + '.lengths'
        if not os.path.exists(args.lens):
            print(f"Lengths file {args.lens} not found, and no lengths file provided", file=sys.stderr)
            sys.exit(1)
        
    return args

def find_chr(starts, lengths):
    offsets = np.cumsum(lengths)
    ### label each mum with the contig it belongs to
    contig_idx = np.searchsorted(offsets, starts, side='right')
    ### get the relative offset of each mum to the start of its contig
    left_start = np.hstack((0, offsets[:-1]))
    rel_offsets = starts - left_start[contig_idx]
    return contig_idx, rel_offsets


def get_contig_names(lengths_file):
    ### assumes lengths_file is formatted as multilengths
    names = []
    cur_name = []
    first_line = True
    for l in open(lengths_file, 'r').readlines():
        l = l.strip().split()
        if first_line and l[1] != '*':
            print('Lengths file must be formatted as multilengths.', file=sys.stderr)
            sys.exit(1)
        first_line = False
        if l[1] == '*':
            if cur_name:
                names.append(cur_name)
            cur_name = []
            continue
        cur_name.append(l[1])
    names.append(cur_name)
    return names

def plot(genome_lengths, polygons, colors, centering, xlims = None, size=None, genomes=None):
    fig, ax = plt.subplots()
    max_length = max(genome_lengths)
    # Just plot simple genome lines
    for idx, g in enumerate(genome_lengths):
        ax.plot([centering[idx] + 0, centering[idx] + g], [idx, idx], 
                alpha=0.2, linewidth=0.75, c='black')
     
    if xlims is not None:
        ax.set_xlim(*xlims) 
    else:
        ax.set_xlim(0, max_length)
    ax.add_collection(PolyCollection(polygons, linewidths=0, alpha=0.8, edgecolors=colors, facecolors=colors))
    
    ax.yaxis.set_ticks(list(range(len(genome_lengths))))
    ax.tick_params(axis='y', which='both',length=0)
    if genomes:
        ax.set_yticklabels(genomes)
    else:
        ax.yaxis.set_ticklabels([])
    
    ax.set_xlabel('genomic position')
    ax.set_ylabel('sequences')
    ax.set_ylim(-0.25, len(genome_lengths)-1 + 0.25)
    # ax.invert_yaxis()
    fig.set_tight_layout(True)
    if size:
        fig.set_size_inches(*size)
    return fig, ax


def plot_extract(args, coords, mums, mum_bounds, other_coords, seq_idx, sequences, seq_lengths):
    offsets = np.array(other_coords)
    plot_mums = MUMdata.from_arrays(
        mums.lengths[mum_bounds[0]:mum_bounds[1] + 1].copy(), 
        mums.starts[mum_bounds[0]:mum_bounds[1] + 1, sequences].copy(), 
        mums.strands[mum_bounds[0]:mum_bounds[1] + 1, sequences].copy()
    )
    plot_mums.starts -= offsets[:,0]
    start, end = np.array(coords) - offsets[seq_idx,0]
    centering= [0] * len(sequences)
    poly, colors = mumemto.viz_mums.get_mum_polygons(plot_mums, centering, inv_color='green')
    fig, ax = plot(np.array(seq_lengths)[sequences], poly, colors, centering, xlims=(0, max([x[1] - x[0] for x in other_coords])), size=(10,5))
    ax.plot([start, start], [seq_idx - 0.5, seq_idx + 0.5], color='red', linestyle='--', linewidth=1)
    ax.plot([end, end], [seq_idx - 0.5, seq_idx + 0.5], color='red', linestyle='--', linewidth=1)
    
    # print('saving plot to', os.path.join(args.output, 'extract_synteny.pdf'), file=sys.stderr)
    fig.savefig(args.output + '_extract_synteny.pdf')

def plot_full_synteny(args, coords, mums, mum_bounds, other_coords, seq_idx, sequences, seq_lengths):
    offsets = np.array(other_coords)
    plot_mums = MUMdata.from_arrays(
        mums.lengths.copy(), 
        mums.starts[:, sequences].copy(), 
        mums.strands[:, sequences].copy(),
        blocks = mums.blocks.copy()
    )
    # plot_mums.starts -= offsets[:,0]
    start, end = np.array(coords)# - offsets[seq_idx,0]
    centering= [0] * len(sequences)
    poly, colors = mumemto.viz_mums.get_block_polygons(plot_mums.blocks, plot_mums, centering, inv_color='green')
    fig, ax = plot(np.array(seq_lengths)[sequences], poly, colors, centering, size=(10,5))
    ax.plot([start, start], [seq_idx - 0.5, seq_idx + 0.5], color='red', linestyle='--', linewidth=1)
    ax.plot([end, end], [seq_idx - 0.5, seq_idx + 0.5], color='red', linestyle='--', linewidth=1)
    for i, (start_coord, end_coord) in enumerate(other_coords):
        ax.plot([start_coord, start_coord], [i - 0.5, i + 0.5], color='black', linestyle=':', linewidth=1)
        ax.plot([end_coord, end_coord], [i - 0.5, i + 0.5], color='black', linestyle=':', linewidth=1)
    # print('saving plot to', os.path.join(args.output, 'extract_synteny.pdf'), file=sys.stderr)
    fig.savefig(args.output +'_full_synteny.pdf')

def convert_local_to_global_coords(coords, names, lengths):
    ## coords of format contig:start-end
    coords = coords.split(':')
    contig = coords[0]
    start, end = int(coords[1].split('-')[0]), int(coords[1].split('-')[1])
    assert contig in names, f"sequence {contig} not found in indicated FASTA file"
    offset = sum(lengths[:names.index(contig)])
    return offset + start, offset + end
def convert_global_to_local_coords(start, end, names, lengths):
    contig, rel_offsets = find_chr((start, end), lengths)
    assert contig[0] == contig[1], f"start and end coords are in different contigs: {names[contig[0]]} and {names[contig[1]]}"
    return names[contig[0]], rel_offsets
    
def find_target_region(mums, coords, seq_idx, sequences):    
    collinear_mums = [i for s, e in mums.blocks for i in range(s, e+1)]
    order = np.argsort(mums.starts[collinear_mums, seq_idx])
    coll_mums = MUMdata.from_arrays(mums.lengths[collinear_mums][order], mums.starts[collinear_mums][order], mums.strands[collinear_mums][order])
    starts = coll_mums.starts
    left_mum_idx = bisect_left(starts[:, seq_idx], coords[0]) - 1
    right_mum_idx = bisect_left(starts[:, seq_idx] + coll_mums.lengths, coords[1])
    mum_bounds = (left_mum_idx, right_mum_idx)
    ### deal with case where the requested region is within a mum on either side
    left_mum, right_mum = coll_mums[mum_bounds[0]], coll_mums[mum_bounds[1]]
    left_bound = left_mum.starts[seq_idx]
    right_bound = right_mum.starts[seq_idx]
    left_offset, right_offset = 0, 0
    if coords[0] < left_mum.starts[seq_idx] + left_mum.length:
        left_offset = coords[0] - left_mum.starts[seq_idx]
    if coords[1] > right_mum.starts[seq_idx]:
        right_offset = coords[1] - right_mum.starts[seq_idx]
    print('left margin:', coords[0] - left_bound - left_offset, file=sys.stderr)
    print('right margin:', right_bound + right_offset - coords[1], file=sys.stderr)
    other_coords = [(starts[mum_bounds[0], i] + left_offset, starts[mum_bounds[1], i] + right_offset) for i in sequences]
    return coll_mums, mum_bounds, other_coords

def extract_fasta(args, contig_names, seq_lengths_multi, other_coords, sequences):
    paths = get_seq_paths(args.lens)
    for i in range(len(sequences)):
        p = paths[i]
        with open(p, 'r') as f:
            seq = ''.join(line.strip() for line in f if not line.startswith('>'))
            name, rel_offsets = convert_global_to_local_coords(other_coords[i][0], other_coords[i][1], contig_names[i], seq_lengths_multi[i])
            coord_line = f"{name}:{rel_offsets[0]}-{rel_offsets[1]}"
            with open(os.path.join(args.output, os.path.basename(p).replace('.fa', f'.extract.fa')), 'w') as out:
                out.write(f'>{os.path.splitext(os.path.basename(p))[0]}_{coord_line}\n{seq[other_coords[i][0] : other_coords[i][1]]}\n')

def extract_bed(args, contig_names, seq_lengths_multi, other_coords, sequences):
    paths = get_seq_paths(args.lens)
    bed_file = open(args.output + '.bed', 'w')
    for i in range(len(sequences)):
        p = paths[i]
        name, rel_offsets = convert_global_to_local_coords(other_coords[i][0], other_coords[i][1], contig_names[i], seq_lengths_multi[i])
        bed_file.write(f'{name}\t{rel_offsets[0]}\t{rel_offsets[1]}\t{p}\n')

def main(args=None):
    args = parse_arguments(args)
    if args.fasta and not os.path.exists(args.output):
        os.makedirs(args.output)
    seq_lengths_multi = get_sequence_lengths(args.lens, multilengths=True)
    seq_lengths = [sum(x) for x in seq_lengths_multi]
    contig_names = get_contig_names(args.lens)
    NUM_SEQS = len(seq_lengths)
    if args.sequences is not None and any([s >= NUM_SEQS for s in args.sequences]):
        print(f"Sequence index {max(args.sequences)} is invalid (N = {NUM_SEQS})", file=sys.stderr)
        sys.exit(1)
    if args.sequences is None:
        args.sequences = list(range(NUM_SEQS))
        
    mums = MUMdata(args.mum_file)
    if mums.blocks is None:
        mums.blocks = find_coll_blocks(mums)
        print('found blocks:', len(mums.blocks), file=sys.stderr)
    else:
        print('using provided blocks:', len(mums.blocks), file=sys.stderr)
    coords = convert_local_to_global_coords(args.range, contig_names[args.seq_idx], seq_lengths_multi[args.seq_idx])
    coll_mums, mum_bounds, other_coords = find_target_region(mums, coords, args.seq_idx, args.sequences)
    if args.plot:
        if not args.lens:
            raise ValueError("--lengths-file is required when --plot is specified")
        plot_extract(args, coords, coll_mums, mum_bounds, other_coords, args.seq_idx, args.sequences, seq_lengths)
    if args.plot_full:
        plot_full_synteny(args, coords, mums, mum_bounds, other_coords, args.seq_idx, args.sequences, seq_lengths)
    if args.fasta:
        extract_fasta(args, contig_names, seq_lengths_multi, other_coords, args.sequences)
    else:
        extract_bed(args, contig_names, seq_lengths_multi, other_coords, args.sequences)
    

if __name__ == "__main__":
    main()
