#!/bin/bash

# extract sequence of promoters

REFSEQDIR='data/refseqv1.0/161010_Chinese_Spring_v1.0_pseudomolecules.fasta'
SIMPLEBED='data/refseqv1.1/refseqv1.1_promoters_all_module_genes.bed'
#SIMPLEBED='data/refseqv1.1/refseqv1.1_promoters_two_one_module_genes.bed'
OUTFASTA='data/refseqv1.1/refseqv1.1_promoters_all_module_genes.fa'
#OUTFASTA='data/refseqv1.1/refseqv1.1_promoters_two_one_module_genes.fa'

bedtools getfasta -name -s -fi $REFSEQDIR -bed $SIMPLEBED -fo $OUTFASTA