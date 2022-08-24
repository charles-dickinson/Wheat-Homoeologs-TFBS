#!/bin/bash

# run FIMO

PROMOTERSEQ='data/refseqv1.1/refseqv1.1_promoters_all_module_genes.fa'
#PROMOTERSEQ='data/refseqv1.1/refseqv1.1_promoters_two_one_module_genes.fa'

BACKGROUND='data/refseqv1.1/refseqv1.1_promoters_all_module_genes.fa_background'
#BACKGROUND='data/refseqv1.1/refseqv1.1_promoters_two_one_module_genes.fa_background'

MOTIF_FILE='data/plantPAN/PWM_all.meme'
FIMO_OUTPUT='data/fimo_results/fimo_v1.1_all_module_genes.tsv'
#FIMO_OUTPUT='data/fimo_results/fimo_v1.1_two_one_module_genes.tsv'


fasta-get-markov $PROMOTERSEQ $BACKGROUND

fimo --bgfile $BACKGROUND --skip-matched-sequence --text --no-qvalue $MOTIF_FILE $PROMOTERSEQ > $FIMO_OUTPUT