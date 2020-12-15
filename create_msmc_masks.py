#!/usr/bin/env python3


""" Create bed file with covered regions for demographic inference

Coverage is considered insufficient if it is smaller than the mean of the
coverage core distribution minus 1.5 times the standard deviation of the core
distribution; or higher than the mean plus 1.5 times the standard deviation.

"""

import argparse
import gzip
import io
import pysam

from Bio import SeqIO


def get_args():

    parser = argparse.ArgumentParser(description='Create bed files indicating \
                                     called and mapable regions in the genome.')

    parser.add_argument("-ref",
                        required=True,
                        help='Reference genome in fasta.')

    parser.add_argument('-bam',
                        required=True,
                        help='Bam file for which the masks shall be generated')

    parser.add_argument('-chr', nargs= '+',
                        help='Chromosome name if masks are only required \
                        for a single chromosome')

    parser.add_argument('-cov_m', required=True,
                        help='Median coverage', type=int)

    parser.add_argument('-cov_s', required=True,
                        help='Standard deviation of coverage', type=int)

    parser.add_argument('-baseq',
                        default=20, type=int,
                        help='Minimum base quality.')

    parser.add_argument('-mapq',
                        default=20, type=int,
                        help='Minimum mapping quality.')

    args = parser.parse_args()
    return args


class MaskGenerator:
    """ From msmc-tools (utils.py)
    """

    def __init__(self, filename, chr):
        self.lastCalledPos = -1
        self.lastStartPos = -1
        self.file = io.TextIOWrapper(gzip.open(filename, "w"))
        self.chr = chr

    # assume 1-based coordinate, output in bed format
    def addCalledPosition(self, pos):
        if self.lastCalledPos == -1:
            self.lastCalledPos = pos
            self.lastStartPos = pos
        elif pos == self.lastCalledPos + 1:
            self.lastCalledPos = pos
        else:
            self.file.write("{}\t{}\t{}\n".format(self.chr, self.lastStartPos - 1, self.lastCalledPos))
            self.lastStartPos = pos
            self.lastCalledPos = pos


def main():

    args = get_args()
    ref = args.ref
    baseq, mapq = args.baseq, args.mapq
    bam = args.bam

    cov_median, cov_stdev = args.cov_m, args.cov_s

    contigs = [(seq_record.id, len(seq_record))
               for seq_record in SeqIO.parse(ref, "fasta")]

    acc_name = bam.split('/')[-1].split('.')[0]

    #bad_positions = {}


    for c in contigs:

        chromosome, start, end = c[0], 0, c[1]
        if args.chr and chromosome not in args.chr:
            continue

        print("Creating masks on chromosome %s ..." % (chromosome))

        #bad_positions[chromosome] = {}

        outname = "%s.%s.mask.bed.gz" % (acc_name, chromosome)
        mask = MaskGenerator(outname, chromosome)
        pybam = pysam.AlignmentFile(bam, "rb")

        for pileupcolumn in pybam.pileup(chromosome, start, end, **{"truncate": True}):

            pos = int(pileupcolumn.pos + 1)
            n_good, n_bad = 0, 0

            for pileupread in pileupcolumn.pileups:

                read = pileupread.alignment

                if pileupread.is_del or pileupread.is_refskip:
                    continue
                elif (read.mapping_quality < mapq) or (read.query_qualities[pileupread.query_position] < baseq) or not read.is_paired:
                    n_bad += 1
                else:
                    n_good += 1


            #outline = [chromosome, str(pos), str(n_good), str(n_bad)]
            #out.write('\t'.join(outline) + '\n')

            tot_cov = n_good + n_bad
            if n_good > 0 and (cov_median-1.5*cov_stdev) <= tot_cov <= (cov_median+1.5*cov_stdev) and (n_bad / float(n_good)) < 1:
                mask.addCalledPosition(pos)

            #else:
            #    bad_positions[chromosome][pos] = (n_good, n_bad)

        pybam.close()


if __name__ == '__main__':
    main()
