#!/usr/bin/env python3

import argparse
import sys, os
import numpy as np
from mumemto.utils import MUMdata, get_sequence_lengths, get_seq_paths
import mumemto.viz_mums
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
from tqdm.auto import tqdm

def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description="Shred a MUM file into smaller fragments and optionally visualize.")
    parser.add_argument("mum_file", type=str, help="Path to the input MUM file.")
    parser.add_argument("--shred-size", type=int, default=10_000_000, help="Size of each shred in base pairs (default: 10Mbp).")
    parser.add_argument("--plot", action="store_true", help="If set, generate a visualization of the shreds.")
    parser.add_argument("--fasta", action="store_true", help="If set, generate a FASTA file for each shred. Otherwise, only write a BED file.")
    parser.add_argument("--output", '-o', type=str, default="shreds", help="Output directory for shredded MUMs.")
    parser.add_argument("--sequences", '-s', type=int, nargs='*', default=None, help="One or more sequence indices to output BED or FASTA for. By default, all sequences are included.")
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

def plot(genome_lengths, polygons, colors, centering, size=None, genomes=None, filename=None):
    fig, ax = plt.subplots()
    max_length = max(genome_lengths)
    # Just plot simple genome lines
    for idx, g in enumerate(genome_lengths):
        ax.plot([centering[idx] + 0, centering[idx] + g], [idx, idx], 
                alpha=0.2, linewidth=0.75, c='black')
      
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
    ax.set_xlim(0, max_length)
    fig.set_tight_layout(True)
    ax.invert_yaxis()

    if size:
        fig.set_size_inches(*size)
    return fig, ax

def shred_mums(mums, shred_size):
    NUM_SEQS = mums.strands.shape[1]
    shred_starts = np.zeros((NUM_SEQS, 1), dtype=int)
    
    for _, e in mums.blocks:
        if np.all(mums.strands[e,:]) and (mums.starts[e,:] + mums.lengths[e] - shred_starts[:,-1]).mean() > shred_size:
            shred_starts = np.column_stack((shred_starts, mums.starts[e,:] + mums.lengths[e]))
    
    return mums, shred_starts

def shred_fasta(args, shredded_mums):
    NUM_SHREDS = shredded_mums.shape[1]
    paths = get_seq_paths(args.lens)
    output_files = [open(os.path.join(args.output, f'shred{i}.fa'), 'w') for i in range(NUM_SHREDS)]
    for idx in args.sequences:
        p = paths[idx]
        with open(p, 'r') as f:
            seq = ''.join(line.strip() for line in f if not line.startswith('>'))
            # assert len(seq) == seq_lengths[idx], f'length mismatch: {p}, expected {seq_lengths[idx]}, got {len(seq)}'
            for i in range(NUM_SHREDS - 1):
                output_files[i].write(f'>{p}\n{seq[shredded_mums[idx, i] : shredded_mums[idx, i+1]]}\n')
            output_files[-1].write(f'>{p}\n{seq[shredded_mums[idx, -1] : ]}\n')
    for f in output_files:
        f.close()
      
def shred_bed(args, shredded_mums):
    NUM_SHREDS = shredded_mums.shape[1]
    seq_lengths = get_sequence_lengths(args.lens)
    names = get_contig_names(args.lens)
    for i in args.sequences:
        contig_idx, rel_offsets = find_chr(shredded_mums[i,:], seq_lengths[i])
        with open(os.path.join(args.output, f'sequence{i}.bed'), 'w') as f:
            for j in range(NUM_SHREDS - 1):
                f.write(f'{names[contig_idx[j]]}\t{rel_offsets[j]}\t{rel_offsets[j + 1]}\n')
def plot_shreds(args, mums, shred_starts):
    NUM_SEQS = mums.strands.shape[1]
    seq_lengths = get_sequence_lengths(args.lens)
    centering = [0] * NUM_SEQS
    
    poly, colors = mumemto.viz_mums.get_block_polygons(mums.blocks, mums, centering, inv_color='green')
    fig, ax = plot(seq_lengths, poly, colors, centering, size=(10,15))
    
    # Plot the shred boundaries
    for i in range(shred_starts.shape[1]):
        ax.plot(shred_starts[:,i], np.arange(NUM_SEQS), color='red', linewidth=2)
    
    fig.savefig(os.path.join(args.output, 'shred_synteny.pdf'))

def main(args=None):
    args = parse_arguments(args)
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    seq_lengths = get_sequence_lengths(args.lens)
    NUM_SEQS = len(seq_lengths)
    if args.sequences is not None and any([s >= NUM_SEQS for s in args.sequences]):
        print(f"Sequence index {max(args.sequences)} is invalid (N = {NUM_SEQS})", file=sys.stderr)
        sys.exit(1)
    if args.sequences is None:
        args.sequences = list(range(NUM_SEQS))
    mums = MUMdata(args.mum_file)
    mums, shredded_mums = shred_mums(mums, args.shred_size)
    if args.fasta:
        shred_fasta(args, shredded_mums)
    else:
        shred_bed(args, shredded_mums)
    if args.plot:
        if not args.lens:
            raise ValueError("--lengths-file is required when --plot is specified")
        plot_shreds(args, mums, shredded_mums)

if __name__ == "__main__":
    main()
