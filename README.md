Tools for population genetics
==============================

## detectIEH.py
Detect islands of extended heterozygosity (IEH) among called SNPs, as they arise
when reads from paralogous stretches of DNA are collapsed onto a single region
in the reference genome. Such "fake" heterozygosity often survives conventional
SNP hard-filtering and can bias downstream analyses.

## pairwise_nucdiff.py
Calculate pairwise genetic differences (dXY) from a vcf file which contains
variant and non-variant sites.

## n-fold_degenerate.py
Given a reference genome sequence and a gene annotation in gff format, output
positions with a codon degeneracy of n (4 by default).
