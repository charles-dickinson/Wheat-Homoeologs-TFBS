#!/bin/bash

# make bed file containing positions of 1.5kb upstream of TSS

REFSEQGFF=data/refseqv1.0/refseqv1.0.gff3
CHRENDPOINTS=data/refseqv1.0/chromosome_endpoints.txt
BEDDIR=data/refseqv1.0/refseqv1.0.bed
OUTDIR=data/refseqv1.0/refseqv1.0_promoters.bed

gff2bed < $REFSEQGFF > $BEDDIR
bedtools flank -i $BEDDIR -g $CHRENDPOINTS -l 1500 -r 0 -s > $OUTDIR
