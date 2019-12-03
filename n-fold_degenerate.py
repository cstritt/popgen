#!/usr/bin/env python2
# -*- coding: utf-8 -*-


""" Extract fourfold degenerate site positions from annotated reference genome

Input:
    - reference genome
    - gff annotation with cds phases (1-based index!)

"""

import re
import argparse
import sys

from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio import SeqIO
from Bio.Seq import Seq


def get_args():

    parser = argparse.ArgumentParser(description='Get positions of\
                                     n-fold degenerate sites in an annotated\
                                     reference genome.')
    parser.add_argument("-ref",
                        dest="reference",
                        required=True,
                        help='Reference genome in fasta format.')

    parser.add_argument('-gff',
                        dest="annot",
                        required=True,
                        help='Gene annotation in gff format. Has to contain \
                        codon positions.')

    parser.add_argument('-n',
                        dest="degeneracy",
                        type = int, default = 4,
                        help='Codon degeneracy. 4D sites by default.')

    args = parser.parse_args()
    return args


class annotation_entry:
    def __init__(self, line):

        fields = line.strip().split('\t')

        self.chromosome = fields[0]
        self.start = int(fields[3])
        self.end = int(fields[4])
        self.id = re.search(r'ID=(.*?)[;\n]', line).group(1)
        self.feature = fields[2]
        self.strand = fields[6]
        self.phase = fields[7]


def altcodons(codon, table):
    """List codons that code for the same aminonacid / are also stop.

    @param codon
    @table code table id
    @return list of codons

    """
    tab = unambiguous_dna_by_id[table]

    if codon in tab.stop_codons:
        return tab.stop_codons
    try:
        aa = tab.forward_table[codon]
    except:
        return []

    return [k for (k, v) in tab.forward_table.iteritems()
            if v == aa and k[0] == codon[0] and k[1] == codon[1]]


def load_files(reference, annot):

    ref = {}
    for seq_record in SeqIO.parse(reference, "fasta"):
        ref[seq_record.id] = seq_record.seq

    d = {}
    with open(annot) as f:
        for line in f:
            if line.startswith('#'):
                continue
            annot = annotation_entry(line)
            if annot.feature == 'CDS':

                parent = '.'.join(annot.id.split('.')[:-2])
                if annot.chromosome not in d:
                    d[annot.chromosome] = {}

                if parent not in d[annot.chromosome]:
                    d[annot.chromosome][parent] = []
                d[annot.chromosome][parent].append(annot)
    return ref, d


def full_cds(gene, ref):

    seq_concat = Seq('')
    for cds in gene:
        if cds.strand == '-':
            start = cds.start - 1
            end = cds.end
            seq = ref[cds.chromosome][start:end].reverse_complement()
        else:
            start = cds.start
            seq = ref[cds.chromosome][start:end]

        seq_concat += seq
    return seq_concat


def fourfold_degenerate_positions(codons):

    out = []

    for c in codons:

        alt = altcodons(str(c), 1)

        first = len({x[0] for x in alt})
        second = len({x[1] for x in alt})
        third = len({x[2] for x in alt})

        out += [first, second, third]

    return out


def get_degenerate_positions(ref, d, n):

    out = {'chromosome':  [],
           'position': [],
           'gene': []}

    for chrmsm in d:
        for gene in d[chrmsm]:

            seq = Seq('')
            posizioni = []  # zero-based!
            strand = d[chrmsm][gene][0].strand

            for cds in d[chrmsm][gene]:

                start = cds.start - 1
                end = cds.end

                if cds.strand == '-':
                    posizioni_cds = range(start, end)[::-1]
                    seq_cds = ref[chrmsm][start:end].reverse_complement()
                else:
                    posizioni_cds = range(start, end)
                    seq_cds = ref[chrmsm][start:end]

                seq += seq_cds
                posizioni += posizioni_cds

            codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
            fdp = fourfold_degenerate_positions(codons)

            """ Convert 0-based to 1-based by adding 1 if the gene is on the
            forward strand or leaving it if on the reverse
            """
            if strand == '-':
                degen = [posizioni[i] for i, x in enumerate(fdp) if x == n]
            else:
                degen = [posizioni[i]+1 for i, x in enumerate(fdp) if x == n]

            out['chromosome'] += len(degen) * [chrmsm]
            out['position'] += degen
            # out['gene'] += len(degen) * [gene]

    tabula = [[out['chromosome'][i], str(out['position'][i])]
              for i in range(len(out['chromosome']))]

    tabula.sort(key=lambda x: (x[0], int(x[1])))

    return tabula


def main():

    args = get_args()
    ref, d = load_files(args.reference, args.annot)
    tabula = get_degenerate_positions(ref, d, args.degeneracy)
    for line in tabula:
        sys.stdout.write('\t'.join(line) + '\n')


if __name__ == '__main__':
    main()
