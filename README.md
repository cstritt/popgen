Tools for population genetics
==============================

### Methods used in the paper 'Migration without interbreeding: evolutionary history of a highly selfing grass inferred from whole genomes'
C Stritt, EL Gimmi, M Wyler, AH Bakali, A Skalska, R Hasterok, LAJ Mur, N Pecchioni, AC Roulin

https://www.biorxiv.org/content/10.1101/2020.09.03.280842v2

#### detectIEH.py
Detect islands of extended heterozygosity (IEH) among called SNPs, as they arise
when reads from paralogous stretches of DNA are collapsed onto a single region
in the reference genome. Such "fake" heterozygosity often survives conventional
SNP hard-filtering and can bias downstream analyses.

#### pairwise_nucdiff.py
Calculate pairwise genetic differences (dXY) from a vcf file which contains
variant and non-variant sites.

#### n-fold_degenerate.py
Given a reference genome sequence and a gene annotation in gff format, output
positions with a codon degeneracy of n (4 by default).

#### create_msmc_masks.py
Create bed file indicating regions within 1.5 standard deviations of the mean
of the coverage. Used for demographic inference with MSMC2 and to obtain high
confidence regions in the genome.
