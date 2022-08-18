#!/bin/bash

# run FIMO

PROMOTERSEQ='data/refseqv1.0/refseqv1.0_promoters.fa'
BACKGROUND='data/refseqv1.0/refseqv1.0_promoters.fa_background'
MOTIF_FILE='data/plantPAN/PWM_all.meme'
FIMO_OUTPUT='data/fimo_results/fimo_v1.0.tsv'

fasta-get-markov $PROMOTERSEQ $BACKGROUND

fimo --bgfile $BACKGROUND --skip-matched-sequence --text --no-qvalue $MOTIF_FILE $PROMOTERSEQ > $FIMO_OUTPUT
