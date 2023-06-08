#!/usr/bin/env python

import sys
import argparse
import pysam, gzip
from Bio import SeqIO
from itertools import chain
from multiprocessing import Pool

class HelpFormatter(
        argparse.RawDescriptionHelpFormatter,
        argparse.ArgumentDefaultsHelpFormatter):
    pass

def read_fa(path):
    handle = gzip.open(path, "rt") if path.endswith(('gz', 'gzip', 'GZ', 'GZIP')) else open(path, 'r')
    for record in SeqIO.parse(handle, "fasta"):
        yield (record.id, str(record.seq).upper())

def out_high_depth_regions(depths, name, seq, min_depth, min_length):
    s = e = -1
    total_valid = 0
    for i, d in enumerate(chain(*depths)):
        if d >= min_depth:
            if s == -1:
                s = e = i
            else:
                e = i
        elif s > -1:
            if e - s + 1 >= min_length:
                print(">%s_%d_%d\n%s" % (name, s, e, seq[s : e + 1]))
                total_valid += e - s + 1
            s = e = -1
    if s > -1:
        if e - s + 1 >= min_length:
            print(">%s_%d_%d\n%s" % (name, s, e + 1, seq[s : e + 1]))
            total_valid += e - s + 1
    print("output rate in %s: %.3f%%" % (name, 100 * total_valid/len(seq)), file=sys.stderr)

def worker(args):
    name, start, end, bamfile = args #end is not include
    depth = [0] * (end - start)
    bam_reader = pysam.AlignmentFile(bamfile)
    for read in bam_reader.fetch(name, start, end):
        if read.is_unmapped or (read.query_alignment_end - read.query_alignment_start) / read.infer_read_length() < 0.8:
            continue
        for i in range(read.reference_start, read.reference_end):
            if start <= i < end:
                depth[i - start] += 1
    return depth

def main(args):
    for name, seq in read_fa(args.genome):
        pool = Pool(args.thread)
        batch_len = int(len(seq)/args.thread) + 1
        batch_len_regions = []
        for depth in pool.imap(worker, [(name, batch_len * i, batch_len * (i + 1), args.bam) for i in range(args.thread)]):
            batch_len_regions.append(depth)
        out_high_depth_regions(batch_len_regions, name, seq, args.min_depth, args.min_len)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class = HelpFormatter,
        description = ''' 
remove_low_depth_in_fasta:
    output sequences with high mapping depth based on a bam file.
    indels in bam does not affect filtering, because the mapping 
    depth is accumulated by the regions of alignments spanned.

exmples: 
    %(prog)s sgs.map.sort.bam genome.fa > genome.filter.fa

'''
    )
    parser.add_argument('-v', '--version', action='version', version='0.1')
    parser.add_argument('bam', 
        help='read-to-ref mapping file in sorted BAM format, index is required.')
    parser.add_argument('genome',
        help='genome assembly file in [GZIP] FASTA format.')
    parser.add_argument('-t', '--thread', metavar = 'INT', type=int, default=5,
        help='number of threads.')
    parser.add_argument('-d', '--min_depth', metavar = 'INT', type=int, default=3,
        help='minimum mapping depth to output.')
    parser.add_argument('-l', '--min_len', metavar = 'INT', type=int, default=1000,
        help='minimum length to output.')
    args = parser.parse_args()
    main(args)
