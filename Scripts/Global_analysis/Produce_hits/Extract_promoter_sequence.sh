#!/bin/bash

# extract sequence of promoters

REFSEQDIR='data/refseqv1.0/161010_Chinese_Spring_v1.0_pseudomolecules.fasta'
SIMPLEBED='data/refseqv1.0/refseqv1.0_promoters_simple.bed'
OUTFASTA='data/refseqv1.0/refseqv1.0_promoters.fa'

bedtools getfasta -name -s -fi $REFSEQDIR -bed $SIMPLEBED -fo $OUTFASTA
