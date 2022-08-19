#!/bin/bash

# get only unique motif hits

FIMO_RESULTS='data/fimo_results/fimo_v1.0.tsv'
FIMO_UNIQUE='data/fimo_results/fimo_v1.0_unique.tsv'

cut -f 1,2,3 $FIMO_RESULTS | uniq > $FIMO_UNIQUE
