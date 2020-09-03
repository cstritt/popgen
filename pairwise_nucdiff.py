#!/usr/bin/env python
# -*- coding: utf-8 -*-


""" Calculate pairwise genetic differences from vcf which contains variant and non-variant sites

Created on Wed Oct 31 11:21:04 2018
@author: cristobal
"""

import argparse
import gzip
import copy
import sys


def get_args():

    parser = argparse.ArgumentParser(description='Count number of pairwise genetic differences')

    parser.add_argument('-v', dest='vcf_file', required=True)

    parser.add_argument('-p', dest='popfile', required=True, help='List of samples to be compared.')

    parser.add_argument('-d', dest='covfile', help='Tab separated file with sampleID - coverage mean -  coverage stdev')

    parser.add_argument('-chr', required=True, nargs = '+', help='Consider specific chromosome/s')

    parser.add_argument('-ignore_het', action='store_true', help='Ignore heterozygous sites')

    parser.add_argument('-m', dest='masks', help='bed file indicating regions which should be considered')

    args=parser.parse_args()
    return args


def pairwise_comparisons(*pops):

    """ Create list with pairwise comparisons. If 2 populations are
    provided, only comparisons between them are returned.
    """
    if len(pops) == 1:
        accessions = pops[0]
        pw = [(i,j) for i in accessions for j in accessions if i != j]
    else:
        accessions_a, accessions_b = pops
        pw = [(i,j) for i in accessions_a for j in accessions_b]
    # remove duplicate pairs (i,j) (j,i)
    for entry in pw:
        if (entry[1], entry[0]) in pw:
            pw.remove((entry[1],entry[0]))
    # set counter for each comparison: (differences, missing data)
    pw_d = {k: [0, 0] for k in pw}
    return pw_d


def create_mask_dictionary(bed_file, chromosomes):
    """ Cave: in bed files, the start coordinate is 0-based,
    the end coordinate 1-based (in contrast to 1-based gff)
    """

    d = {}
    with gzip.open(bed_file) as f:
        next(f) # skip header
        for line in f:
            line = line.decode('utf-8')

            chrm, strt, end = line.strip().split('\t')
            if chrm not in chromosomes:
                continue
            if chrm not in d:
                d[chrm] = set()
            for i in range(int(strt), int(end)+1):
                d[chrm].add(int(i))
    for chrm in d:
        sys.stderr.write('Considering ' + str(len(d[chrm])) + ' positions on ' + chrm + '\n')

    return d

#%%
def main():

    # arguments
    args = get_args()
    vcf_file = args.vcf_file
    popfile = args.popfile
    covfile = args.covfile

    ignore_het = args.ignore_het
    chromosomes = args.chr

    if args.masks:
        masks = create_mask_dictionary(args.masks, chromosomes)
    else:
        masks = []

    # Hard filters. Not advised unless samples all have similar sequencing depths
    #filters = {'GQ' : 0, 'RGQ' : 0, 'min_depth' : 5, 'max_depth' : 50}

    stats = {'multiallelic_or_indel' : 0,
             'bad_qual' : 0,
             'bad_depth' : 0,
             'missing_gt' : 0,
             'het_gt' : 0,
             'masked' : 0}

    with open(popfile) as f:
        pop = [line.strip() for line in f]

    covd = {}
    with open(covfile) as f:
        for line in f:
            acc, cov, cov_std = line.strip().split()
            covd[acc] = (int(cov), int(cov_std))

    chrmsm_d = {}
    counter = int

    with gzip.open(vcf_file) as f:
        for line in f:
            line = line.decode('utf-8')

            if line.startswith('#CHROM'):
                    accs = line.strip().split('\t')[9:]
                    accs_pop = [accs.index(x) for x in accs if x in pop]
                    pop_pwc = pairwise_comparisons(accs_pop)

            elif line.startswith('#'):
                continue

            else:
                fields = line.strip().split('\t')

                chrmsm, pos = fields[:2]

                if chrmsm not in chromosomes:
                    continue

                if not chrmsm in chrmsm_d:
                    counter = 1e6
                    chrmsm_d[chrmsm] = copy.deepcopy(pop_pwc)
                    sys.stderr.write('Working on chromosome ' + chrmsm + '\n')

                if masks and int(pos) not in masks[chrmsm]:
                    stats['masked'] += 1
                    continue

                if int(pos) > counter:
                    sys.stderr.write(pos + '\n')
                    counter += 1e6

                alt = fields[4].split(',')
                if len(alt) > 1 and alt[0] not in ['A', 'C', 'T', 'G']:
                    stats['multiallelic_or_indel'] += 1
                    continue

                info = fields[8].split(':')
                genotypes = [x.split(':') for x in fields[9:]]

                geno_d = []
                for i, gt in enumerate(genotypes):

                    d = {}

                    if len(gt) != len(info):
                        d['GT'] = './.'
                        geno_d.append(d)

                    else:
                        for j in range(len(info)):
                            d[info[j]] = gt[j]
                        geno_d.append(d)


                for comp in chrmsm_d[chrmsm]:

                    a = geno_d[comp[0]]
                    b = geno_d[comp[1]]

                    # Missing genotypes?
                    if a['GT'] == './.' or b['GT'] == './.':
                        stats['missing_gt'] += 1
                        continue

                    a_gt_set = set(a['GT'].split('/'))
                    b_gt_set = set(b['GT'].split('/'))


                   # Heterozygous genotypes?
                    if ignore_het and (len(a_gt_set) > 1 or len(b_gt_set) > 1):
                        stats['het_gt'] += 1
                        continue

                    # FILTERING
                    """ Filtering based on sequencing depth: only count sites which in both a and b
                    have an allele count within 1.5 standard deviations of the mean coverage
                    """

                    a_cov = covd[accs[comp[0]]]
                    b_cov = covd[accs[comp[1]]]

                    if (a_cov[0] - 1.5 * a_cov[1]) < int(a['DP']) < (a_cov[0] + 1.5 * a_cov[1]) and \
                        (b_cov[0] - 1.5 * b_cov[1]) < int(b['DP']) < (b_cov[0] + 1.5 * b_cov[1]):

                        diff = len(a_gt_set - b_gt_set)
                        chrmsm_d[chrmsm][comp][0] += diff
                        chrmsm_d[chrmsm][comp][1] += 1

                    else:
                        stats['bad_depth'] += 1


    for s in stats:
        sys.stderr.write(s + ':' + str(stats[s]) + '\n')

    out = [['chr', 'a', 'b', 'S', 'total_sites', 'dxy', ]]

    for chrmsm in chrmsm_d:
        for comp in chrmsm_d[chrmsm]:

            ind1 = accs[comp[0]]
            ind2 = accs[comp[1]]

            ndiff = chrmsm_d[chrmsm][comp][0]
            ncalled = chrmsm_d[chrmsm][comp][1]

            if ncalled == 0:
                dxy = 'NA'
            else:
                dxy = ndiff / float(ncalled)

            out.append([chrmsm, ind1, ind2, str(ndiff), str(ncalled), str(dxy)])


    for line in out:
        sys.stdout.write('\t'.join(line) + '\n')


if __name__ == '__main__':
    main()
