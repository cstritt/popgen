#!/usr/bin/env python
# -*- coding: utf-8 -*-


""" Detect islands of fake heterozygosity due to copy number variants

Input: vcf file
Output:
    - bed file for each accession in the vcf, indicating suspicuous regions
    - summary for each accession: Nr of het sites, nr of IFHs, mean size of IFHs

"""

import argparse
import gzip
import re
import sys


#%% CLUSTERING


def get_args():

    parser = argparse.ArgumentParser(description='Identify islands of fake heterozygosity')

    parser.add_argument('-f', '--vcf',
                        required=True,
                        help='Gzipped vcf file.')

    parser.add_argument('-c', '--chr',
                        required=False, nargs='+',
                        help='Only consider specific chromosome/s.')

    parser.add_argument('-a', '--acc',
                        required=False, nargs='+',
                        help='Only consider specific accession.')


    args=parser.parse_args()
    return args


class cluster:

    """ Identify clusters of heterozygous sites
    """

    def __init__(self):
        self.sites = []
        self.alt_counts = []
        self.alt_qual = []
        self.ref_counts = []
        self.ref_qual = []
        self.gt_likelihoods = []

    def add_site(self, site):
        self.sites.append(site[0])
        self.ref_counts.append(int(site[1]))
        self.alt_counts.append(int(site[2]))
        #self.alt_qual.append(int(info['QA']))
        #self.ref_qual.append(int(info['QR']))
        #self.gt_likelihoods.append([float(x) for x in info['GL'].split(',')])

    def summarize(self):
        n_sites = len(self.sites)
        start, end = min(self.sites), max(self.sites)
        mean_alt_cov = sum(self.alt_counts) / float(len(self.alt_counts))
        mean_ref_cov = sum(self.ref_counts) / float(len(self.ref_counts))
        return [n_sites, start, end, round(mean_alt_cov), round(mean_ref_cov)]



def create_het_dict(vcf, chromosomes):

    """ Create a dictionary containing the heterozyous sites for each individual, or for a
    subset of accessions and chromosomes
    """

    contigs = []

    with gzip.open(vcf, 'r') as f:

        for line in f:

            # read chromosome lengths
            line = line.decode('utf-8')
            if line.startswith('##contig'):
                name = re.search(r'ID=(.*?),', line).group(1)
                contigs.append(name)


            # read accession names
            elif line.startswith('#CHROM'):
                accs = line.strip().split('\t')[9:]

                het_d = {}

                for a in accs:
                    het_d[a] = {c: [] for c in contigs}

            elif line.startswith('#'):
                continue

            else:
                fields = line.strip().split('\t')
                chromosome = fields[0]
                if chromosomes and chromosome not in chromosomes:
                    continue

                # only biallelic loci
                ref, alt = fields[3:5]
                if len(alt) > 1:
                    continue

                het_d = add_het_site(fields, accs, het_d)

    return het_d, contigs


def prop_missing(genotypes):
    missing_count = 0
    for gt in genotypes:
        if '.' in gt:
            missing_count += 1
    missing_prop = missing_count / float(len(genotypes))
    return missing_prop


def add_het_site(fields, accs, het_d):

    """ Detect clusters of heterozygous sites, indicating
    fake heterozygosity due to reference deletions or non-reference
    insertions
    """

    chromosome = fields[0]
    pos = fields[1]

    format_ids = fields[8].split(':')
    genotypes_raw = fields[9:]

    for i, g in enumerate(genotypes_raw):

        format_l = genotypes_raw[i].split(':')
        format_d = {format_ids[j]: format_l[j] for j in range(len(format_l))}

        acc = accs[i]
        #if acc not in selection:
        #    continue

        gt = set(format_d['GT'].split('/'))

        if len(gt) != 2:
            continue

        if 'AD' in format_d:
            # Some heterozygous sites don't have the AD annotation. GATK bug?
            allele_depth = format_d['AD'].split(',')
            if len(allele_depth) == 2:

                ref_cov, alt_cov = allele_depth
                het_d[acc][chromosome].append([int(pos), ref_cov, alt_cov])

    return het_d


def find_ifhs(het_d, chrmsm, acc,
              max_dist = 200,
              min_n_sites = 10,
              min_length = 300):

    """ Find islands of fake heterozygosity (IFHs).

    max_dist: maximum distance between heterozygous SNPs in an island
    min_n_sites: minimum number of het SNPs in an island
    min_length: minimum length of an island

    """

    resultati = []

    previous = None
    switch = 0

    for site in het_d[acc][chrmsm]:

        pos = site[0]

        if not previous:
            previous = pos
            continue

        if pos - previous < max_dist:

            # new cluster
            if switch == 0:
                cl = cluster()
                switch = 1

            cl.add_site(site)

        else:
            # end cluster
            if switch == 1:
                n_sites, start, end, mean_alt_cov, mean_ref_cov = cl.summarize()

                if end-start > min_length and n_sites > min_n_sites:

                    resultati.append([chrmsm, start, end, (end-start+1), n_sites, mean_ref_cov, mean_alt_cov])

                switch = 0

        previous = pos

    return resultati


def write_output(resultati, acc, chrm):
    """ Write bed file for acc and summary statistics of the bed
    entries.

    Summary: nr of islands, nr of het SNPs in islands,
    min island size, max island size, mean island size,
    min nr het SNPs per island, max nr het SNPs, mean nr het SNPs
    """


    # Write bed
    with open('%s.%s.fake_het.bed' % (acc, chrm), 'w') as f:
        header = ['chrom', 'start', 'end', 'length', 'n_het_sites', 'mean_alt_cov', 'mean_ref_cov']
        f.write('\t'.join(header) + '\n')

        for line in resultati:

            line = list(map(str, line))
            f.write('\t'.join(line) + '\n')


    # Calulcate summary statistics
    n_islands = len(resultati)
    n_sites = 0
    sites_per_island = []
    island_sizes = []

    if not resultati:
        outline = [acc, chrm] + 8*[0]

    else:

        for line in resultati:

            n_sites += line[4]
            sites_per_island.append(line[5])
            island_sizes.append(line[3])

        outline = [acc,
                   chrm,
                   n_islands,
                   n_sites,

                   min(island_sizes),
                   max(island_sizes),
                   (sum(island_sizes) / len(island_sizes)),

                   min(sites_per_island),
                   max(sites_per_island),
                   (sum(sites_per_island) / len(sites_per_island))]

    return outline


#%% MAIN

def main():

    args = get_args()

    vcf = args.vcf
    chromosomes = args.chr
    accessions = args.acc

    het_d, contigs = create_het_dict(vcf, chromosomes)
    if not accessions:
        accessions = het_d.keys()

    summary = [['acc',
                'chr',
                'nr_ifh',
                'nr_het_snps',
                'ifh_min_len',
                'ifh_max_len',
                'ifh_mean_len',
                'nr_snps_min',
                'n_snps_max',
                'n_snps_mean']]

    for acc in accessions:

        sys.stderr.write('Writing output for %s\n' % (acc))

        for c in contigs:

            if chromosomes and c not in chromosomes:
                continue

            res = find_ifhs(het_d, c, acc)
            smry = write_output(res, acc, c)
            summary.append(smry)

    with open('fakeHet_summary.tsv', 'w') as f:
        for line in summary:
            line = list(map(str, line))
            f.write('\t'.join(line) + '\n')


if __name__ == '__main__':
    main()
