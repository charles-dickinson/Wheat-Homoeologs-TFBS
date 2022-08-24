#!/bin/bash

# getting promoters FOR ALL WHEAT GENES

REFSEQGFF=data/refseqv1.1/refseqv1.1.gff3
CHRENDPOINTS=data/refseqv1.1/chromosome_endpoints.txt
BEDDIR=data/refseqv1.1/refseqv1.1.bed
OUTDIR=data/refseqv1.1/refseqv1.1_promoters.bed

gff2bed < $REFSEQGFF > $BEDDIR
bedtools flank -i $BEDDIR -g $CHRENDPOINTS -l 1500 -r 0 -s > $OUTDIR